function [A_phi_sn, A_phi_sp, b_phi_sn, b_phi_sp] = ...
phi_s_mats(p, jn, jp, I, phi_sn_ref, phi_sp_ref)
% phi_s_mats formulates matrices for charge in solid phase.
% The charge in electrolyte is subject to the molar flux J. 
% Note that reference charge at points N- and Nsep are required
% to make the governing equation well defined.

global Nn Np

delta_xn = p.L_n/(Nn-1);
delta_xp = p.L_p/(Np-1);

const_n = p.sig_n*p.epsilon_s_n/delta_xn^2;
const_p = p.sig_p*p.epsilon_s_p/delta_xp^2;

%% formulate A_phi_sn
main_diag_n = -2*ones(1,Nn-1);
upper_diag_n = [2, ones(1,Nn-3)];
lower_diag_n = [ones(1,Nn-3), 2];

A_phi_sn = diag(main_diag_n)+diag(upper_diag_n,1)+diag(lower_diag_n,-1);

A_phi_sn = const_n*[A_phi_sn; [zeros(1,Nn-2),1]];

%% formulate A_phi_sp
main_diag_p = -2*ones(1,Np-1);
upper_diag_p = [2, ones(1,Np-3)];
lower_diag_p = [ones(1,Np-3), 2];

A_phi_sp = diag(main_diag_p)+diag(upper_diag_p,1)+diag(lower_diag_p,-1);

A_phi_sp = const_p*[[1,zeros(1,Np-2)]; A_phi_sp];

%% formulate b_phi_sn
bn_1 = p.a_s_n*p.Faraday*jn(1,1)-2*I/delta_xn;
bn_mid = p.a_s_n*p.Faraday*jn(1,2:end)';
bn_end = const_n*phi_sn_ref;

b_phi_sn = [bn_1; bn_mid; bn_end];

%% formulate b_phi_sp
bp_1 = const_p*phi_sp_ref;
bp_mid = p.a_s_p*p.Faraday*jp(1,1:end-1)';
bp_end = p.a_s_p*p.Faraday*jp(1,end)-2*I/delta_xp;

b_phi_sp = [bp_1; bp_mid; bp_end];

end

