function jac_mat = remove_mass_matrices(jac_mat)
% cremove_mass_matrices remvoes mass matrices on the diagonals of the
% jacobian matrix

global Ns Nn Np Nsep Nccn Nccp

N_csn = Ns*Nn;

% set mass matrices of cs to zero
for i = 1:Nn
    jac_mat((i-1)*Ns+1:i*Ns,(i-1)*Ns+1:i*Ns) = zeros(Ns,Ns);
end

for i = 1:Np
    jac_mat(N_csn+(i-1)*Ns+1:N_csn+i*Ns,N_csn+(i-1)*Ns+1:N_csn+i*Ns) = zeros(Ns,Ns);
end

% set mass matrices of ce, phi_e, phi_s, and T to zero
% Note that df_jdj may be dependent on SEI resistance and thus is not set to zero 
range = Ns*(Nn+Np)+1:Ns*(Nn+Np)+2*(Nn+Nsep+Np)+Nn+Np+Nccn+Nn+Nsep+Np+Nccp;

for i = 1:length(range)
    position = range(i);
    jac_mat(position, position-1:position+1) = 0;
end

end