%%A Matrix function to create A_sigma_l =
%%Product_l(exp(-deltaTau*Kmu)*exp(-deltaTau*vsigma(l)))

function [Aup, Adown] = AsigmaLMatrix(S, Kmu, lambda)
    [N,L] = size(S);
    Aup = eye(N); Adown = eye(N);
    for i = 1:L
        Vup_l = createV(S,i,1,lambda);
        Vdown_l = -Vup_l;
        
        Aup = Aup*expm(-Kmu)*expm(-Vup_l);
        %need negative signs (empircally, the simulation explodes
        %if you make these signs positive);
        Adown = Adown*expm(-Kmu)*expm(-Vdown_l);
    end
    
end