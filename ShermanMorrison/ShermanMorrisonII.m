%% WE HAVE TO BE CAREFUL...the Gubernatis book has exp(V)exp(K) so 
%% all the formula orientations are messed up...we've fixed it below.


function [Gup, Gdown] = ShermanMorrisonII(Gup, Gdown, S, i, l,S0, lambda)
   [N,L] = size(S);
    Iup = eye(N) - Gup; 
    Idown = eye(N) - Gdown;
    change = S(i,l) - S0(i,l);
    
    %%deltaup and delta down are the only places which encode new info
    %%about HS field
    deltaup = zeros(N,N);
    deltadown = zeros(N,N);
    deltaup(i,i) = exp(-2*lambda*S0(i,l))-1; %l indexes deltaTau, the time step we are at
    deltadown(i,i) = exp(2*lambda*S0(i,l))-1;
            
    R_up = 1 + deltaup(i,i)*(1 - Gup(i,i));
    R_down = 1 + deltadown(i,i)*(1 - Gdown(i,i));
    
    %Gchanges = [];
  
    deltaGup = Iup* deltaup * (Gup) / (R_up);
    deltaGdown = Idown*deltadown*(Gdown) / (R_down);
    %Gchanges = [Gchanges; deltaGup, deltaGdown];

    Gup = Gup - deltaGup;
    Gdown = Gdown - deltaGdown;

end