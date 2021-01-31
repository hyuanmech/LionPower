function jac_FD = jacobian_finiteDifference(X, func)

% central difference: eps^(1/3); forward difference: eps^(1/2)
epsilon = eps^(1/3);

M_P = repmat(X, 1, length(X)) + diag(epsilon*ones(length(X),1));
M_N = repmat(X, 1, length(X)) - diag(epsilon*ones(length(X),1));

F_P = splitapply(func, M_P, 1:size(M_P,2));
F_N = splitapply(func, M_N, 1:size(M_N,2));

jac_FD = (F_P-F_N)/2/epsilon;

end