% The Green's matrix, which we derive from:
%Z = SUM(det(Md)det(Mup)), Gup and Gdown are Md and Mup

function [Gup,Gdown] = GreenMatrixUDR(deltaTau,Kmu,S, lambda)
    
    [N,L] = size(S);
    
    v_1up = createV(S,1,lambda);
    B1 = expm(-deltaTau*Kmu)*expm(v_1up); %put kmu as argument
    [Ql,R1, P1] = qr(B1); 
    Dl = diag(diag(R1));
    Tl = Dl^-1*R1*P1.';
    Tproduct = Tl;
    for l = 2:L

       v_lup = createV(S,l,lambda);
       B_l = expm(-deltaTau*Kmu)*expm(v_lup);
       C_l = (B_l*Ql)*Dl; 
       [Ql,R,P] = qr(C_l); %q*r*p.' = C...the original matrix, got to transpose P
       Dl = diag(diag(R));
       Tl = Dl^-1 * R * P.';
       Tproduct = Tl*Tproduct;

    end

    Db = DecompDb(Dl);
    Ds = DecompDs(Dl);
    Gup = (Db^-1*Ql.' + Ds*Tproduct)^-1*(Db^-1*Ql.');
    
    %% Now do the same for Gdown
    v_ldown = -createV(S,1,lambda);
    B1d = expm(-deltaTau*Kmu)*expm(v_ldown);

    [Ql,R1, P1] = qr(B1d); 
    Dl = diag(diag(R1));
    Tl = Dl^-1*R1*P1.';
    Tproduct = Tl;
    for l = 2:L

       v_ldown = -1*createV(S,l,lambda);
       B_l = expm(-deltaTau*Kmu)*expm(v_ldown);
       C_l = (B_l*Ql)*Dl; 
       [Ql,R,P] = qr(C_l); %q*r*p.' = C...the original matrix, got to transpose P
       Dl = diag(diag(R));
       Tl = Dl^-1 * R * P.';
       Tproduct = Tl*Tproduct;

    end
    Db = DecompDb(Dl);
    Ds = DecompDs(Dl);
    Gdown = (Db^-1*Ql.' + Ds*Tproduct)^-1*(Db^-1*Ql.');
    
    
    
    
end