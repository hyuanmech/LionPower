function D_s_p = Dsp_calc(p, T)
% compute diffusion coefficient in particles of positive electrode

global Ns

D_s_p = 1.25e-15*exp(-p.E_D_p/p.R*(1/T-1/298.15));
D_s_p = D_s_p*ones(Ns,1);

end