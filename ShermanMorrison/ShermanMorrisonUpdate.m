function [Gup, Gdown, c_kup, b_jup] = ...
    ShermanMorrisonUpdate(Gup, Gdown, S, i, l,S0, lambda)

   [N,L] = size(S);

    up = exp(-2*lambda*S0(i,l))-1;
    down = exp(2*lambda*S0(i,l))-1;
    
    %%what index ordering should we do for Gup(j,i) or Gup(i,j)
    %indexing of G should be (G(i,j)*G(k,i)
    c_kup = -(up)*Gup(:,i);
    c_kup(i) = c_kup(i) + up;
    c_kdown =  -(down)*Gdown(:,i);
    c_kdown(i) = c_kdown(i) + down;

    b_jup = Gup(i,:)/(1 + c_kup(i));
    b_jdown = Gdown(i,:)/(1 + c_kdown(i));

    Gup = Gup - c_kup*b_jup;
    Gdown = Gdown - c_kdown*b_jdown;

end