function k_p = kp_calc(p, T)
% compute rate constant of positive electrode
% based on [Li et al., Power Source 255 (2014) 130-143]

k_p = 1.4e-12*exp(-p.E_k_p/p.R*(1./T-1/298.15));

end