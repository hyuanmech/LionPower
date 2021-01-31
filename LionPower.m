%==============================================================%
% LionPower: An Open-Source Model for Electrochemical-Thermal  % 
% Analysis of Lithium-ion Battery                              % 
%                                                              % 
% Copyright (c) 2021 Hao Yuan                                  % 
%==============================================================%

clc;
clear;
close all;

global Ns Nn Nsep Np Nt Nccn Nccp Tref I dt_all

%% Model setup
Ns = 11;
Nn = 11;
Nsep = 4;
Nccn = 10;
Nccp = 10;
Np = Nn;
Tref = 298.15;
flag_solver = 0; % 0: default fsolve; 1: jacobian free Newton Krylov
if flag_solver == 0
    fprintf('<strong>LionPower fsolve</strong>\n');
else
    fprintf('<strong>LionPower JFNK</strong>\n');
end

%% Start parallel pool
pc = gcp;

tic;
%% Import parameters for electrochemical model
run Parameters_LFPO_Graphite

%% Experimental data
BatteryCapacity = 12;   % [Ah]
ExpDataName = [{'Discharge_1C'}, {'Discharge_2C'}, {'UDDS'}];
selectData = 2; % choose from 1, 2, and 3
ExpData = xlsread('Measurements.xlsx',ExpDataName{selectData});
t = ExpData(:,1);
current = -ExpData(:,2);
                   
I = current/p.Area;                      % [A/m^2]
Nt = length(t);

dt_all = t - [0; t(1:end-1)];
t_out = [0; t];

%% Initial condition and variables for storage
SOC = 100; % SOC at the beginning
Tini = 298.15;
[InitialVar, StateVar] = BatteryInitialization(p, SOC, Tini);

csn0 = InitialVar{1,1};
csp0 = InitialVar{2,1};
rn = InitialVar{3,1};
rp = InitialVar{4,1};
phi_e0 = InitialVar{5,1};

csn_all = StateVar{1,1};
csp_all = StateVar{2,1};
cs_surf_n_all = StateVar{3,1};
cs_surf_p_all = StateVar{4,1};
ce_all = StateVar{5,1};
phi_e_all = StateVar{6,1};
phi_sn_all = StateVar{7,1};
phi_sp_all = StateVar{8,1};
jn_all = StateVar{9,1};
jp_all = StateVar{10,1};
T_all = StateVar{11,1};
q_tot = StateVar{12,1};
V_out = StateVar{13,1};

%% Compute phi_e, phi_s, and j when current is firstly applied
% Note that ce and cs do not change at the very beginning.
fprintf('Computing consistent initial conditions:\n');
flag = 0;
if mean(I) ~= 0
    I_mean = abs(mean(I));
else
    I_mean = abs(mean(abs(max(I)), abs(min(I))));
end
jn_mean = I_mean/p.a_s_n/p.Faraday/p.L_n;
jp_mean = I_mean/p.a_s_p/p.Faraday/p.L_p;

% Mean value of j is used as initial guess for fast convergence 
jn_all(:,1) = sign(I(1))*jn_mean*ones(Nn,1);
jp_all(:,1) = -sign(I(1))*jp_mean*ones(Np,1);

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
fprintf('\tFormulating the initial Jacobian matrix ...');
n_tot_init = (Nn+Nsep+Np)+2*(Nn+Np);
x = sym('X', [n_tot_init,1]);
F_residual_init = BatteryEqnsAll(x, args_init, flag);
jac_mat_init = jacobian(F_residual_init, x);

% build jacobian function
jacFuncName_init = 'jac_mat_init';
jac_init = jacobian_postprocessing(jac_mat_init, jacFuncName_init);
fprintf(' Done!\n');

% solve the equations at time 0
fprintf('\tSolving the equations ...');
options = optimoptions(@fsolve, 'Display', 'none',...
    'SpecifyObjectiveGradient', true, 'CheckGradients', false,...
    'Diagnostics','off', 'FunctionTolerance',1e-12,...
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
fprintf(' Done!\n');

%% Start simulation
fprintf('Main simulation: \n');
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

if flag_solver == 0
    % formulate the jacobian matrix
    fprintf('\tFormulating the full Jacobian matrix ...');
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
    fprintf(' Done!\n');
end

fprintf('\tSolving the equations at t = %.4f s', t(1));
for hh = 2:Nt
        
    % The values at the previous step are taken as the initial guess
    X = [csn_all(:,hh-1)/p.c_s_n_max;...
        csp_all(:,hh-1)/p.c_s_p_max;...
        ce_all(:,hh-1)/p.c_e;...
        phi_e_all(:,hh-1);...
        phi_sn_all(:,hh-1);...
        phi_sp_all(:,hh-1);...
        T_all(:,hh-1)/Tref;...
        jn_all(:,hh-1)/jn_mean;...
        jp_all(:,hh-1)/jp_mean;...
        I(hh-1)/I_mean];

    % parameters passed to discretised PDE equations
    args = cell(7,1);
    args{1,1} = p;
    args{2,1} = T_all(:,hh-1)/Tref;
    args{3,1} = csn_all(:,hh-1)/p.c_s_n_max;
    args{4,1} = csp_all(:,hh-1)/p.c_s_p_max;
    args{5,1} = ce_all(:,hh-1)/p.c_e;
    args{6,1} = I(hh)/I_mean;
    args{7,1} = hh;
    
    if flag_solver == 0
        
        % solve the equations
        f = @(X) BatteryEqnsAll(X, args, flag, jac);
        
        % default fsolve
        options = optimoptions(@fsolve, 'Display', 'none',...
        'SpecifyObjectiveGradient', true, 'CheckGradients', false,...
        'Diagnostics','off', 'FunctionTolerance',1e-6,...
        'FiniteDifferenceType','central', 'FunValCheck','on',...
        'UseParallel',false);

        [X,F] = fsolve(f,X,options);  
        
        [~, jac_all{hh-1,1}] = BatteryEqnsAll(X, args, flag, jac);
        
    else
        
        % Jacobian Free Newton Krylov Method 
        % solve the equations
        f = @(X) BatteryEqnsAll(X, args, flag);
        
        if hh == 2
            % construct preconditioner for Newton-Krylov method
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
        end
        X_prev = X;

    end
    
    N_csn = Ns*Nn;
    N_csp = Ns*Np;
    N_ce = Nn+Nsep+Np;
    N_phi_e = Nn+Nsep+Np;
    N_phi_sn = Nn;
    N_phi_sp = Np;
    N_T = Nccn+Nn+Nsep+Np+Nccp;
    N_jn = Nn;
    N_jp = Np;

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
    
    strTime_prev = num2str(t(hh-1), '%.4f');
    strLen = length(strTime_prev);
    strTime = num2str(t(hh), '%.4f');
    fprintf([repmat('\b',1,strLen+2), strTime, ' s']);
        
    if phi_sp_all(end,hh)-phi_sn_all(1,hh)<2.5
        t_out(hh+1:end) = [];
        break;
    end

end
fprintf(' Done!\n');

save(strcat('LionPower_', ExpDataName{selectData}, '_T', ...
    num2str(p.T_amb-273.15), '_', num2str(BatteryCapacity/p.Area),'.mat'), ...
    't_out', 'csn_all', 'csp_all', 'ce_all', 'phi_e_all', ...
    'phi_sn_all', 'phi_sp_all', 'T_all', 'jn_all', 'jp_all');

toc;