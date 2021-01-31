function [U,dUrefdT] = PotentialAnode(theta, T)
% Potential for negative electrode: Unref(theta_n, T)
% Obtained from [Safari and Delacourt, J. Electrochem. Soc. 2011 A562-A571]
% The temperature-dependent potential is approximated by Taylor's first
% order expansion

Uref = 0.6379+0.5416*exp(-305.5309*theta)...
    +0.044*tanh(-(theta-0.1958)/0.1088)...
    -0.1978*tanh((theta-1.0571)/0.0854)...
    -0.6875*tanh((theta+0.0117)/0.0529)...
    -0.0175*tanh((theta-0.5692)/0.0875);

dUrefdT = (344.1347148*exp(-32.9633287*theta+8.316711484)...
    ./(1+749.0756003*exp(-34.79099646*theta+8.887143624))...
    -0.8520278805*theta+0.362299229*theta.^2+0.2698001697)/1000;

U = Uref + (T-298.15).*dUrefdT;
    
end

