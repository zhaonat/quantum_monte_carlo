%% function to create K matrix for the PAM
%% there is no filling potential for the uncorrelated electrons...

function K = createK_l(Nx, Ny, Nf, V_k, t, mu)
    HoppingPart = KE2DMatrix(Nx,Ny, t);
    ChemicalPotential = mu*eye(Nf*Nf);
    
    K = blkdiag(HoppingPart, ChemicalPotential);
    
    %% Now we incorporate all the interaction terms
    

end
