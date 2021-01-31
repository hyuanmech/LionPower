function [U,dUrefdT] = PotentialCathode(theta, T)
% Potential for positive electrode: Upref(theta_p, T)
% Obtained from [Safari and Delacourt, J. Electrochem. Soc. 2011 A562-A571]
% The temperature-dependent potential is approximated by Taylor's first
% order expansion

Uref = 3.4323-0.8428*exp(-80.2493*(1-theta).^1.3198)...
    -3.2474e-6*exp(20.2645*(1-theta).^3.8003)...
    +3.2482e-6*exp(20.2646*(1-theta).^3.7995);

dUrefdT = -0.35376*theta.^8+1.3902*theta.^7-2.2585*theta.^6 ...
    +1.9635*theta.^5-0.98716*theta.^4+0.28857*theta.^3 ...
    -0.046272*theta.^2+0.0032158*theta-1.9186e-5;

U = Uref + (T-298.15).*dUrefdT;

    
end