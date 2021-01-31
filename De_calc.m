function D_e = De_calc(c_e, T)
% compute diffusion coefficient of Li salt in electrolyte

c_e = c_e/1000;

D_e = 10.^(-8.43-54./(T-229-5*c_e)-0.22*c_e);

% D_e = 1.5e-10*ones(size(c_e));

% D_e = 10.^(-4.43-(54./(Tmean-229-0.005*c_e))-0.22*0.001*c_e-4);


end