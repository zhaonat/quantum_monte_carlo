% The Green's matrix, which we derive from:

%Z = SUM(det(Md)det(Mup)), Gup and Gdown are Md and Mup
function [Gup, Gdown] = GreenMatrixLPerm(A,S, lambda,Lperm)
    [N,L] = size(S);
    Gup = eye(N); Gdown = eye(N);
    
    for l = Lperm
        Vup_l = createV(S,l,lambda);
        Vdown_l = -Vup_l;
        
        Gup = Gup*A*expm(Vup_l);
        Gdown = Gdown*A*expm(Vdown_l);
    end
    
    %% add identity 
    Gup = eye(N) + Gup;
    Gdown = eye(N) + Gdown;
    
    %% invert
    Gup = Gup^-1;
    Gdown = Gdown^-1;
    
end