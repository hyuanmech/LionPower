function kappa = kappa_calc(ce, T)
% compute conductivity of electrolyte

% [Jiang et al., PS, 2013]
% c_e = c_e/1000; % convert [mol/m3] to [mol/L]
% kappa = 20.8409*c_e-21.29*c_e.^2+13.6986*c_e.^3-7.58544*c_e.^4+...
%     2.45464*c_e.^5-0.30446*c_e.^6;


% % [Northrop et al., ECS, 2011]
% Tmean = mean(T);
% kappa = 1e-4*ce.*((-10.5+0.668*1e-3*ce+0.494*1e-6*ce.^2) +...
%     (0.074  -1.78*1e-5*ce -8.86*1e-10*ce.^2).*Tmean + ...
%     (-6.96*1e-5+2.8*1e-8*ce).*Tmean.^2).^2;

ce = ce/1000;
kappa = 1.2544*ce.*(-8.2488+0.053248*T-2.9871e-5*T.^2 ...
    +0.26235*ce-9.3063e-3*ce.*T+8.069e-6*ce.*T.^2+0.22002*ce.^2-1.765e-4*ce.^2.*T).^2;

end