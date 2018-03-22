%% index l denotes which time slice 

function V = createV_l(Nx, Ny, S, l, lambda)
    UncorrelatedPart = eye(Nx, Ny);
    CorrelatedPart = lambda*diag(S(:,l));
    V = blkdiag(UncorrelatedPart, CorrelatedPart);
end