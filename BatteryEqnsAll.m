function [F,J] = BatteryEqnsAll(X, args, flag, jacFunc)
% BatteryEqns formulates discretised non-linear equations.
% X is the estimation, which MUST come from previous results to compute b.
% flag MUST be either 0 or 1. 
% flag = 0, initial condition, flag = 1, main simulation

global Nccn Nn Ns Nsep Np Nccp dt_all Tref I

p = args{1,1};
T_prev = args{2,1};
csn_prev = args{3,1};
csp_prev = args{4,1};
ce_prev = args{5,1};
I_scalar = args{6,1};
hh = args{7,1};

delta_t = dt_all(hh-1);

%% ================================== Construct PDE system =====================================
N_csn = Ns*Nn;
N_csp = Ns*Np;
N_ce = Nn+Nsep+Np;
N_phi_e = Nn+Nsep+Np;
N_phi_sn = Nn;
N_phi_sp = Np;
N_T = Nccn+Nn+Nsep+Np+Nccp;
N_jn = Nn;
N_jp = Np;

% cs, ce, jn, and jp have no dimension
if flag == 1
    csn = X(1:N_csn,1);
    csp = X(N_csn+1:N_csn+N_csp,1);
    cs_surf_n = csn(Ns:Ns:end,1);
    cs_surf_p = csp(Ns:Ns:end,1);
    ce = X(N_csn+N_csp+1:N_csn+N_csp+N_ce,1);
    phi_e = X(N_csn+N_csp+N_ce+1:N_csn+N_csp+N_ce+N_phi_e,1);
    phi_sn = X(N_csn+N_csp+N_ce+N_phi_e+1:N_csn+N_csp+N_ce+N_phi_e+N_phi_sn,1);
    phi_sp = X(N_csn+N_csp+N_ce+N_phi_e+N_phi_sn+1:N_csn+N_csp+N_ce+N_phi_e+N_phi_sn+N_phi_sp,1);
    T = X(N_csn+N_csp+N_ce+N_phi_e+N_phi_sn+N_phi_sp+1:...
        N_csn+N_csp+N_ce+N_phi_e+N_phi_sn+N_phi_sp+N_T,1);
    jn = X(N_csn+N_csp+N_ce+N_phi_e+N_phi_sn+N_phi_sp+N_T+1:...
        N_csn+N_csp+N_ce+N_phi_e+N_phi_sn+N_phi_sp+N_T+N_jn,1);
    jp = X(N_csn+N_csp+N_ce+N_phi_e+N_phi_sn+N_phi_sp+N_T+N_jn+1:...
        N_csn+N_csp+N_ce+N_phi_e+N_phi_sn+N_phi_sp+N_T+N_jn+N_jp,1);
    I_density = X(end);
else
    csn = zeros(N_csn,1); % not required for initial condition
    csp = zeros(N_csp,1); % not required for initial condition
    cs_surf_n = csn_prev(Ns:Ns:end,1);
    cs_surf_p = csp_prev(Ns:Ns:end,1);
    ce = zeros(N_ce,1); % not required for initial condition
    phi_e = X(1:N_phi_e,1);
    phi_sn = X(N_phi_e+1:N_phi_e+N_phi_sn,1);
    phi_sp = X(N_phi_e+N_phi_sn+1:N_phi_e+N_phi_sn+N_phi_sp,1);
    T = zeros(N_T,1); % not required for initial condition
    jn = X(N_phi_e+N_phi_sn+N_phi_sp+1:...
        N_phi_e+N_phi_sn+N_phi_sp+N_jn,1);
    jp = X(N_phi_e+N_phi_sn+N_phi_sp+N_jn+1:...
        N_phi_e+N_phi_sn+N_phi_sp+N_jn+N_jp,1);
    I_density = 0; % not required for initial condition
end

% set initial ionic flux to absolute average value
if mean(I) ~= 0
    I_mean = mean(I);
else
    I_mean = mean(abs(max(I)), abs(min(I)))*sign(abs(max(I))-abs(min(I)));
end
jn_mean = I_mean/p.a_s_n/p.Faraday/p.L_n;
jp_mean = I_mean/p.a_s_p/p.Faraday/p.L_p;

if flag == 1
    % equations for csn and csp
    A_csn_all = cell(Nn,1);
    for i = 1:Nn
        [Atn, Fn] = csn_mass_mats_log_grids_CVM_dimensionless(p, T(Nccn+i,1));
        gn = [zeros(Ns-1,1); -jn(i,1)*jn_mean/p.R_s_n/p.c_s_n_max];
        A_csn = 2*Atn-delta_t*Fn;
        b_csn = (2*Atn+delta_t*Fn)*csn_prev((i-1)*Ns+1:i*Ns,1)+2*delta_t*gn;
        f_csn((i-1)*Ns+1:i*Ns,1) = A_csn*csn((i-1)*Ns+1:i*Ns,1)-b_csn;
        A_csn_all{i,1} = A_csn;
    end
    
    A_csp_all = cell(Np,1);
    for i = 1:Np
        [Atp, Fp] = csp_mass_mats_log_grids_CVM_dimensionless(p, T(Nccn+Nn+Nsep+i,1));
        gp = [zeros(Ns-1,1); -jp(i,1)*jp_mean/p.R_s_p/p.c_s_p_max];
        A_csp = 2*Atp-delta_t*Fp;
        b_csp = (2*Atp+delta_t*Fp)*csp_prev((i-1)*Ns+1:i*Ns,1)+2*delta_t*gp;
        f_csp((i-1)*Ns+1:i*Ns,1) = A_csp*csp((i-1)*Ns+1:i*Ns,1)-b_csp;  
        A_csp_all{i,1} = A_csp;
    end
    
    % equations for ce
    [A_ce, b_ce] = ce_mats_FVM_dimensionless(p, ce_prev, jn, jp, jn_mean, jp_mean, T, hh);

    f_ce = A_ce*ce-b_ce;
    
    % equations for T
    q = heat_sources_calc(p, phi_sn, phi_sp, phi_e, ce, ...
        cs_surf_n, cs_surf_p, T, jn, jp, jn_mean, jp_mean, I_density*I_mean);

    [A_T, b_T] = thermal_model_FVM_dimensionless(p, T_prev, q, hh);

    f_T = A_T*T-b_T;
    
    % equation for I_density
    f_I = I_density-I_scalar;
end


% equations for phi_e
[A_phi_e, b_phi_e] = phi_e_mats_FVM_dimensionless(p, ...
    (1-flag)*ce_prev+flag*ce, (1-flag)*T_prev+flag*T, ...
    jn, jp, jn_mean, jp_mean);

f_phi_e = A_phi_e*phi_e-b_phi_e;

% equations for phi_s
[A_phi_sn, A_phi_sp, b_phi_sn, b_phi_sp] = ...
    phi_s_mass_mats_FVM_dimensionless(p, jn, jp, jn_mean, jp_mean,...
    ((1-flag)*I_scalar+flag*I_density)*I_mean);

f_phi_sn = A_phi_sn*phi_sn-b_phi_sn;
f_phi_sp = A_phi_sp*phi_sp-b_phi_sp;

% equations for j
k_n = kn_calc(p, ((1-flag)*T_prev(Nccn+1:Nccn+Np,1)+...
    flag*T(Nccn+1:Nccn+Np,1))*Tref);
k_p = kp_calc(p, ((1-flag)*T_prev(Nccn+Nn+Nsep+1:Nccn+Nn+Nsep+Np,1)+...
    flag*T(Nccn+Nn+Nsep+1:Nccn+Nn+Nsep+Np,1))*Tref);

f_jn = jn-Butler_Volmer_n_dimensionless(phi_sn, p, k_n,...
    (1-flag)*T_prev(Nccn+1:Nccn+Np,1)+flag*T(Nccn+1:Nccn+Np,1),...
    (1-flag)*ce_prev+flag*ce, phi_e, cs_surf_n, jn, jn_mean);
f_jp = jp-Butler_Volmer_p_dimensionless(phi_sp, p, k_p,...
    (1-flag)*T_prev(Nccn+Nn+Nsep+1:Nccn+Nn+Nsep+Np,1)+...
    flag*T(Nccn+Nn+Nsep+1:Nccn+Nn+Nsep+Np,1),...
    (1-flag)*ce_prev+flag*ce, phi_e, cs_surf_p, jp, jp_mean);

if flag == 1
    F = [f_csn; f_csp; f_ce; f_phi_e; f_phi_sn; f_phi_sp; f_T; f_jn; f_jp; f_I];
else
    F = [f_phi_e; f_phi_sn; f_phi_sp; f_jn; f_jp];
end

%% ====================================== JACOBIAN matrix ======================================
if nargout == 2
    
    if flag == 1

        J = zeros(length(X), length(X));

        % mass matrices of csn and csp
        for i = 1:Nn
            J((i-1)*Ns+1:i*Ns,(i-1)*Ns+1:i*Ns) = A_csn_all{i,1};
        end

        for i = 1:Np
            J(N_csn+(i-1)*Ns+1:N_csn+i*Ns,N_csn+(i-1)*Ns+1:N_csn+i*Ns) = A_csp_all{i,1};
        end

        % mass matrix of ce
        J(N_csn+N_csp+1:N_csn+N_csp+N_ce, N_csn+N_csp+1:N_csn+N_csp+N_ce) = A_ce;

        % mass matrix of phi_e
        J(N_csn+N_csp+N_ce+1:N_csn+N_csp+N_ce+N_phi_e, ...
            N_csn+N_csp+N_ce+1:N_csn+N_csp+N_ce+N_phi_e) = A_phi_e;

        % mass matrices of phi_s
        J(N_csn+N_csp+N_ce+N_phi_e+1:N_csn+N_csp+N_ce+N_phi_e+Nn,...
            N_csn+N_csp+N_ce+N_phi_e+1:N_csn+N_csp+N_ce+N_phi_e+Nn) = A_phi_sn;

        J(N_csn+N_csp+N_ce+N_phi_e+Nn+1:N_csn+N_csp+N_ce+N_phi_e+Nn+Np,...
            N_csn+N_csp+N_ce+N_phi_e+Nn+1:N_csn+N_csp+N_ce+N_phi_e+Nn+Np) = A_phi_sp;

        % mass matrix of T
        J(N_csn+N_csp+N_ce+N_phi_e+Nn+Np+1:N_csn+N_csp+N_ce+N_phi_e+Nn+Np+N_T,...
            N_csn+N_csp+N_ce+N_phi_e+Nn+Np+1:N_csn+N_csp+N_ce+N_phi_e+Nn+Np+N_T) = A_T;

        % off-diagonal elements
        if nargin == 4
            J = jacFunc(X, J);
        else
            J(N_csn+N_csp+N_ce+N_phi_e+Nn+Np+N_T+1:end,...
                N_csn+N_csp+N_ce+N_phi_e+Nn+Np+N_T+1:end) = ...
                eye(length(N_csn+N_csp+N_ce+N_phi_e+Nn+Np+N_T+1:length(X)));
        end
    else
        J = zeros(length(X), length(X));
        J = jacFunc(X, J);
    end

end