function [InitialVar, StateVar] = BatteryInitialization(p, SOC, Tini)

global Ns Nn Nsep Np Nt Nccn Nccp

InitialVar = cell(5,1);
StateVar = cell(13,1);

% lithium concentration in electrodes
%------------------------------------------------------
[csn0_value,csp0_value,OCPn,OCPp] = init_cs(p,SOC,Tini);

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

InitialVar{1,1} = csn0;
InitialVar{2,1} = csp0;
InitialVar{3,1} = rn;
InitialVar{4,1} = rp;

StateVar{1,1} = csn_all;
StateVar{2,1} = csp_all;
StateVar{3,1} = cs_surf_n_all;
StateVar{4,1} = cs_surf_p_all;

% lithium concentration in eletrolyte
%------------------------------------------------------
ce0 = p.c_e*ones(Nn+Nsep+Np,1);
ce_all = zeros(Nn+Nsep+Np,Nt);
ce_all(:,1) = ce0;

StateVar{5,1} = ce_all;
%------------------------------------------------------

% charge in electrolyte phase
%------------------------------------------------------
phi_e0 = zeros(Nn+Nsep+Np,1);
phi_e_all = zeros(Nn+Nsep+Np,Nt);
phi_e_all(1,:) = 0; % reference phi_e at 0-

InitialVar{5,1} = phi_e0;
StateVar{6,1} = phi_e_all;
%------------------------------------------------------

% charge in solid phase
%------------------------------------------------------
phi_sn0 = OCPn*ones(Nn,1);
phi_sp0 = OCPp*ones(Np,1);

phi_sn_all = zeros(Nn,Nt);
phi_sp_all = zeros(Np,Nt);

phi_sn_all(:,1) = phi_sn0;
phi_sp_all(:,1) = phi_sp0;

StateVar{7,1} = phi_sn_all;
StateVar{8,1} = phi_sp_all;
%------------------------------------------------------

% molar flux
%------------------------------------------------------
jn_all = zeros(Nn, Nt);
jp_all = zeros(Np, Nt);

StateVar{9,1} = jn_all;
StateVar{10,1} = jp_all;
%------------------------------------------------------

T_all = zeros(Nccn+Nn+Nsep+Np+Nccp,Nt);
T_all(:,1) = Tini*ones(Nccn+Nn+Nsep+Np+Nccp,1);

q_tot = zeros(Nn+Nsep+Np,Nt);

V_out = zeros(1,Nt);

StateVar{11,1} = T_all;
StateVar{12,1} = q_tot;
StateVar{13,1} = V_out;