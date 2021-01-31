function D_s_n = Dsn_calc(p, T)
% compute diffusion coefficient in particles of negative electrode

global Ns

D_s_n = 3.9e-14*exp(-p.E_D_n/p.R*(1/T-1/298.15));
D_s_n = D_s_n*ones(Ns,1);

end