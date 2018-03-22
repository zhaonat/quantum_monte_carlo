% The Green's matrix, which we derive from:

%Z = SUM(det(Md)det(Mup)), Gup and Gdown are Md and Mup
function [Gup, Gdown] = GreenMatrix(deltaTau,Kmu,S, lambda)
    [N,L] = size(S);
    Gup = eye(N); Gdown = eye(N);
    
    for l = 1:L
        Vup_l = createV(S,l,lambda);
        Vdown_l = -Vup_l;
        
        Gup = Gup*expm(-deltaTau*Kmu)*expm(Vup_l); %need negative signs (empircally, the simulation explodes
        %if you make these signs positive);
        %deltaTau not needed for V terms... as it is encoded in lambda
        Gdown = Gdown*expm(-deltaTau*Kmu)*expm(Vdown_l);
    end
    
    %% add identity 
    Gup = eye(N) + Gup;
    Gdown = eye(N) + Gdown;
    
    %% invert
    Gup = Gup^-1;
    Gdown = Gdown^-1;
    
end