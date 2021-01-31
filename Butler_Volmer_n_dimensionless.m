function jn = Butler_Volmer_n_dimensionless(phi_sn, p, k_n, T, ce, phi_e, cs_surf_n, jn, jn_mean)
% Apply Butler-Volmer equation to calculate molar flux j

global Nn Tref

i0n = p.Faraday*k_n.*(ce(1:Nn)*p.c_e).^p.alph...
    .*(p.c_s_n_max-cs_surf_n*p.c_s_n_max).^p.alph.*(cs_surf_n*p.c_s_n_max).^p.alph;

thetan = cs_surf_n*p.c_s_n_max/p.c_s_n_max;

[Un,~] = PotentialAnode(thetan, T*Tref);

etan = phi_sn-phi_e(1:Nn)-Un-p.Faraday*p.R_f_n.*jn*jn_mean;

jn = i0n/p.Faraday.*(exp(p.alph*p.Faraday/p.R./(T*Tref).*etan)-...
    exp(-p.alph*p.Faraday/p.R./(T*Tref).*etan))/jn_mean;

end

