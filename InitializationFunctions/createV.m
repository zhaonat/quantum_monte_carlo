%create potential energy associated with U, the interaction term
%V contains a slice of every column in the HS field (S);
function V = createV(S, l, lambda)

    %why does the matrix to be diagonal (we are basically evaluating a V
    %for every site in the system
    V = lambda*diag(S(:,l));
    
end

