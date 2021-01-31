function jac_mat = jacobian_parallel(F, X)

% number of workers available for the parallel computations
pc = gcp;
n = pc.NumWorkers;

jac_mat = [];

s = floor(length(X)/n);
s_r = length(X)-(n-1)*s;

s_array = [s*ones(1,n-1), s_r];

F_cell = cell(n,1);

% slice F into n sections
for i = 1:n
    F_cell{i,1} = F((i-1)*s+1:(i-1)*s+s_array(i));
end

J_cell = cell(n,1);

% compute the jacobian matrix by sections
parfor i = 1:n
    J_cell{i,1} = jacobian(F_cell{i,1}, X);
end

for i = 1:n
    jac_mat = [jac_mat; J_cell{i,1}];
end

end