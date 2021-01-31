function k_n = kn_calc(p, T)
% compute rate constant of negative electrode
% based on [Li et al., Power Source 255 (2014) 130-143]

k_n = 3e-11*exp(-p.E_k_n/p.R*(1./T-1/298.15));

end