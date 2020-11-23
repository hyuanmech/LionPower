%===========================================
% Electrochemical-thermal model of battery
% Author: Hao Yuan
%===========================================

clc;
clear;
close all;

global dt_all Ns Nn Nsep Np Nt Nccn Nccp Tref I YY

%% Psuedo data for coding
%================================================================
Ns = 11;
Nn = 11;
Nsep = 4;
Nccn = 10;
Nccp = 10;
Np = Nn;
Tref = 298.15;
flag_solver = 0; % 0: default fsolve; 1: jacobian free Newton Krylov

N_csn = Ns*Nn;
N_csp = Ns*Np;
N_ce = Nn+Nsep+Np;
N_phi_e = Nn+Nsep+Np;
N_phi_sn = Nn;
N_phi_sp = Np;
N_T = Nccn+Nn+Nsep+Np+Nccp;
N_jn = Nn;
N_jp = Np;
%=================================================================

%% Import parameters for electrochemical model
run Parameters_LFPO_Graphite_BJTU

%% Experimental data
BatteryCapacity = 12;                    % [Ah]
% load time and current data
load('UDDS_25.mat');
t_array = UDDS_25(:,1);
Current_array = -UDDS_25(:,2); 
t = t_array;                             % [s]
current = Current_array;                 % [A]
                   
I = current/p.Area;                      % [A/m^2]
Nt = length(t);

dt_all = t - [0; t(1:end-1)];

%% Initial condition and variables for storage
% lithium concentration in electrodes
%------------------------------------------------------
SOC = 90; % SOC at the beginning
[csn0_value,csp0_value,OCPn,OCPp] = init_cs(p,SOC);

csn_all = zeros(Ns*Nn, Nt);
csp_all = zeros(Ns*Np, Nt);

cs_surf_n_all = zeros(Nn, Nt);
cs_surf_p_all = zeros(Np, Nt);

csn0 = csn0_value*ones(Ns, Nn-1);
csp0 = csp0_value*ones(Ns, Np-1);

rn = 0:p.R_s_n/(Ns-1):p.R_s_n;
rp = 0:p.R_s_p/(Ns-1):p.R_s_p;

csn_all(:,1) = csn0_value*ones(Ns*Nn,1);
csp_all(:,1) = csp0_value*ones(Ns*Np,1);

cs_surf_n_all(:,1) = csn0_value*ones(Nn,1);
cs_surf_p_all(:,1) = csp0_value*ones(Np,1);


% lithium concentration in eletrolyte
%------------------------------------------------------
ce0 = p.c_e*ones(Nn+Nsep+Np,1);
ce_all = zeros(Nn+Nsep+Np,Nt);
ce_all(:,1) = ce0;
%------------------------------------------------------

% charge in electrolyte phase
%------------------------------------------------------
phi_e0 = zeros(Nn+Nsep+Np,1);
phi_e_all = zeros(Nn+Nsep+Np,Nt);
phi_e_all(1,:) = 0; % reference phi_e at 0-
%------------------------------------------------------

% charge in solid phase
%------------------------------------------------------
phi_sn0 = OCPn*ones(Nn,1);
phi_sp0 = OCPp*ones(Np,1);

phi_sn_all = zeros(Nn,Nt);
phi_sp_all = zeros(Np,Nt);

phi_sn_all(:,1) = phi_sn0;
phi_sp_all(:,1) = phi_sp0;
%------------------------------------------------------

% molar flux
%------------------------------------------------------
jn_all = zeros(Nn, Nt);
jp_all = zeros(Np, Nt);
%------------------------------------------------------

T_all = zeros(Nccn+Nn+Nsep+Np+Nccp,Nt);
Tini = 298.15;
T_all(:,1) = Tini*ones(Nccn+Nn+Nsep+Np+Nccp,1);

q_tot = zeros(Nn+Nsep+Np,Nt);

V_out = zeros(1,Nt);

YY.x = [];
YY.fval = [];

%% Compute phi_e, phi_s, and j when current is firstly applied
% Note that ce and cs do not change at the very beginning.
flag = 0;
if mean(I) ~= 0
    I_mean = abs(mean(I));
else
    I_mean = abs(mean(abs(max(I)), abs(min(I))));
end
jn_mean = I_mean/p.a_s_n/p.Faraday/p.L_n;
jp_mean = I_mean/p.a_s_p/p.Faraday/p.L_p;
 

% The initial guess of phi_e and phi_s come from the static condition
X0 = [phi_e_all(:,1);...
    phi_sn_all(:,1); ...
    phi_sp_all(:,1); ...
    jn_all(:,1)/jn_mean; ...
    jp_all(:,1)/jp_mean];  

% parameters to be passed to discretised PDE equations
args_init = cell(7,1);
args_init{1,1} = p;
args_init{2,1} = T_all(:,1)/Tref;
args_init{3,1} = csn_all(:,1)/p.c_s_n_max;
args_init{4,1} = csp_all(:,1)/p.c_s_p_max;
args_init{5,1} = ce_all(:,1)/p.c_e;
if I(1) == 0 && flag_solver == 1
    nonzero_index = find(I);
    args_init{6,1} = I(nonzero_index(1))/I_mean;
else
    args_init{6,1} = I(1)/I_mean;
end
args_init{7,1} = 2;

% formulate the jacobian matrix at time 0
n_tot_init = (Nn+Nsep+Np)+2*(Nn+Np);
x = sym('X', [n_tot_init,1]);
F_residual_init = BatteryEqnsAll(x, args_init, flag);
jac_mat_init = jacobian(F_residual_init, x);

% build jacobian function
jacFuncName_init = 'jac_mat_init';
jac_init = jacobian_postprocessing(jac_mat_init, jacFuncName_init);

% solve the equations at time 0
options = optimoptions(@fsolve, 'Display', 'iter',...
    'SpecifyObjectiveGradient', true, 'CheckGradients', false,...
    'Diagnostics','on', 'FunctionTolerance',1e-12,...
    'OptimalityTolerance', 1e-12, 'StepTolerance',1e-12,...
    'FiniteDifferenceType','central', 'FunValCheck','on');

f_init = @(X) BatteryEqnsAll(X, args_init, flag, jac_init);

[X0,F0] = fsolve(f_init,X0,options);

% Update phi_e, phi_s, and j
phi_e_all(:,1) = X0(1:Nn+Nsep+Np,1);
phi_sn_all(:,1) = X0(Nn+Nsep+Np+1:Nn+Nsep+Np+Nn,1);
phi_sp_all(:,1) = X0(Nn+Nsep+Np+Nn+1:Nn+Nsep+Np+Nn+Np,1);
jn_all(:,1) = X0(Nn+Nsep+Np+Nn+Np+1:Nn+Nsep+Np+Nn+Np+Nn,1)*jn_mean;
jp_all(:,1) = X0(Nn+Nsep+Np+Nn+Np+Nn+1:Nn+Nsep+Np+Nn+Np+Nn+Np,1)*jp_mean;

%% Start simulation
flag = 1;

% parameters used to build the jacobian matrix
% Since these parameters are updated at each time step,
% the intial values are used for convenience
args = cell(7,1);
args{1,1} = p;
args{2,1} = T_all(:,1)/Tref;
args{3,1} = csn_all(:,1)/p.c_s_n_max;
args{4,1} = csp_all(:,1)/p.c_s_p_max;
args{5,1} = ce_all(:,1)/p.c_e;
args{6,1} = I(1)/I_mean;
args{7,1} = 2;

% formulate the jacobian matrix
n_tot = Ns*(Nn+Np)+3*(Nn+Nsep+Np)+2*(Nn+Np)+Nccn+Nccp+1;
X = sym('X', [n_tot,1]);
F_residual = BatteryEqnsAll(X, args, flag);
jac_mat = jacobian_parallel(F_residual, X);

% % Default jacobian function is slow and thus not used
% jac_mat= jacobian(F_residual, X);

% set all mass matrices in the jacobian matrix to zeros
% Mass matrices are formulated in BatteryEqnsAll and are on the diagonals
% of the jacobian matrix. To save computation time, these matrices are
% inserted into the jacobian matrix directly instead of computing each
% element via the numeric equations.
jac_mat = remove_mass_matrices(jac_mat);

% build jacobian function
jacFuncName = 'jac_mat';
jac = jacobian_postprocessing(jac_mat, jacFuncName);

% debug
jac_all = cell(30,1);

for hh = 2:length(t)
        
    % The values at the previous step are taken as the initial guess
    X = [csn_all(:,hh-1)/p.c_s_n_max;...
        csp_all(:,hh-1)/p.c_s_p_max;...
        ce_all(:,hh-1)/p.c_e;...
        phi_e_all(:,hh-1);...
        phi_sn_all(:,hh-1);...
        phi_sp_all(:,hh-1);...
        T_all(:,hh-1)/Tref;...
        I_mean/I_mean*ones(Nn,1);...
        -I_mean/I_mean*ones(Np,1);...
        I_mean/I_mean];

    % parameters passed to discretised PDE equations
    args = cell(7,1);
    args{1,1} = p;
    args{2,1} = T_all(:,hh-1)/Tref;
    args{3,1} = csn_all(:,hh-1)/p.c_s_n_max;
    args{4,1} = csp_all(:,hh-1)/p.c_s_p_max;
    args{5,1} = ce_all(:,hh-1)/p.c_e;
    args{6,1} = I(hh)/I_mean;
    args{7,1} = hh;
    
    % solve the equations
    f = @(X) BatteryEqnsAll(X, args, flag, jac);
    
    if flag_solver == 0
        
        % default fsolve
        options = optimoptions(@fsolve, 'Display', 'iter',...
        'SpecifyObjectiveGradient', true, 'CheckGradients', false,...
        'Diagnostics','on', 'FunctionTolerance',1e-6,...
        'FiniteDifferenceType','central', 'FunValCheck','on',...
        'UseParallel',false);

        [X,F] = fsolve(f,X,options);  
        
        [~, jac_all{hh-1,1}] = BatteryEqnsAll(X, args, flag, jac);
        
    else
        
        % Jacobian Free Newton Krylov Method 
        if hh == 2
            % construct preconditioner for Newton-Krylov method
            X(end) = I(hh-1)/I_mean;
            X(end-Np:end-1) = jp_all(:,hh-1)/jp_mean;
            X(end-Np-Nn:end-Np-1) = jn_all(:,hh-1)/jn_mean;
            jac_FD = jacobian_finiteDifference(X, f);
            [L,U] = ilu(sparse(jac_FD),struct('type','ilutp','droptol',1e-12));
        end
        
        % GMRES may fail to achieve the specified tolerance, which,
        % however, does not affect the convergence at each time step. 
        warning('off','MATLAB:gmres:tooSmallTolerance');
        
        relTol = 1e-6;
        absTol = 1e-6;
        gmresTol = 1e-3;
        gmresMaxit = 20;
        F = BatteryEqnsAll(X, args, flag);
        r0 = norm(F,2);
        tic;
        X_ini = X;
        gmres_counter = 0; 
        while norm(F,2)>relTol*r0+absTol || sum(isnan(F)) ~= 0
            if isreal(F) == 0 || gmres_counter == 50 || sum(isnan(F)) ~= 0
                X = X_ini;
                X(end) = I(hh)/I_mean;
                X(end-Np:end-1) = -I(hh)/I_mean*ones(Np,1);
                X(end-Np-Nn:end-Np-1) = I(hh)/I_mean*ones(Nn,1);
                F = BatteryEqnsAll(X, args, flag);
                r0 = norm(F,2);
            end
            if gmres_counter == 100
                X = X_prev;
                X(end) = I(hh-1)/I_mean;
                X(end-Np:end-1) = -I(hh-1)/I_mean*ones(Np,1);
                X(end-Np-Nn:end-Np-1) = I(hh-1)/I_mean*ones(Nn,1);
                F = BatteryEqnsAll(X, args, flag);
                jac_FD = jacobian_finiteDifference(X, f);
                [L,U] = ilu(sparse(jac_FD),struct('type','ilutp','droptol',1e-12));
            end
            Jv = @(v) Jv_FD(v, X, F, f);
            [v,fl,relres,iter,resvec] = ...
                gmres(Jv, F, [], gmresTol, gmresMaxit, L, U);        
            X = X - v;
            F = BatteryEqnsAll(X, args, flag); 
            gmres_counter = gmres_counter + 1;
            disp(norm(F));
        end
        toc;
        X_prev = X;
    end
    
    csn_all(:,hh) = X(1:N_csn,1)*p.c_s_n_max;
    csp_all(:,hh) = X(N_csn+1:N_csn+N_csp,1)*p.c_s_p_max;
    cs_surf_n_all(:,hh) = csn_all(Ns:Ns:end,1);
    cs_surf_p_all(:,hh) = csp_all(Ns:Ns:end,1);
    ce_all(:,hh) = X(N_csn+N_csp+1:N_csn+N_csp+N_ce,1)*p.c_e;
    phi_e_all(:,hh) = X(N_csn+N_csp+N_ce+1:N_csn+N_csp+N_ce+N_phi_e,1);
    phi_sn_all(:,hh) = X(N_csn+N_csp+N_ce+N_phi_e+1:N_csn+N_csp+N_ce+N_phi_e+N_phi_sn,1);
    phi_sp_all(:,hh) = X(N_csn+N_csp+N_ce+N_phi_e+N_phi_sn+1:...
        N_csn+N_csp+N_ce+N_phi_e+N_phi_sn+N_phi_sp,1);
    T_all(:,hh) = X(N_csn+N_csp+N_ce+N_phi_e+N_phi_sn+N_phi_sp+1:...
        N_csn+N_csp+N_ce+N_phi_e+N_phi_sn+N_phi_sp+N_T,1)*Tref;
    jn_all(:,hh) = X(N_csn+N_csp+N_ce+N_phi_e+N_phi_sn+N_phi_sp+N_T+1:...
        N_csn+N_csp+N_ce+N_phi_e+N_phi_sn+N_phi_sp+N_T+N_jn,1)*jn_mean;
    jp_all(:,hh) = X(N_csn+N_csp+N_ce+N_phi_e+N_phi_sn+N_phi_sp+N_T+N_jn+1:...
        N_csn+N_csp+N_ce+N_phi_e+N_phi_sn+N_phi_sp+N_T+N_jn+N_jp,1)*jp_mean;
    
    disp(hh);
    
    if phi_sp_all(end,hh)-phi_sn_all(1,hh)<2.5
        break;
    end

end
