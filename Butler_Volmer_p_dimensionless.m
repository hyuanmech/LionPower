function jp = Butler_Volmer_p_dimensionless(phi_sp, p, k_p, T, ce, phi_e, cs_surf_p, jp, jp_mean)
% Apply Butler-Volmer equation to calculate molar flux j

global Nn Nsep Np Tref

i0p = p.Faraday*k_p.*(ce(Nn+Nsep+1:Nn+Nsep+Np)*p.c_e).^p.alph...
    .*(p.c_s_p_max-cs_surf_p*p.c_s_p_max).^p.alph.*(cs_surf_p*p.c_s_p_max).^p.alph;

thetap = cs_surf_p*p.c_s_p_max/p.c_s_p_max;

[Up,~] = PotentialCathode(thetap, T*Tref);

etap = phi_sp-phi_e(Nn+Nsep+1:end)-Up-p.Faraday*p.R_f_p.*jp*jp_mean;

jp = i0p/p.Faraday.*(exp(p.alph*p.Faraday/p.R./(T*Tref).*etap)-...
    exp(-p.alph*p.Faraday/p.R./(T*Tref).*etap))/jp_mean;

end

