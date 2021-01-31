function [csn0,csp0,OCPn,OCPp] = init_cs(p,SOC,Tini)
% Compute Initial Solid Concentrations from Voltage and Params

SOC = SOC/100;      % convert to a fraction
% calculate initial solid concentrations by interpolation
csn0 = ((SOC*(p.theta_max_n-p.theta_min_n) + p.theta_min_n))*p.c_s_n_max;
csp0 = ((SOC*(p.theta_max_p-p.theta_min_p) + p.theta_min_p))*p.c_s_p_max;

theta_n = csn0/p.c_s_n_max;
theta_p = csp0/p.c_s_p_max;

[OCPn,~] = PotentialAnode(theta_n, Tini);
[OCPp,~] = PotentialCathode(theta_p, Tini);
