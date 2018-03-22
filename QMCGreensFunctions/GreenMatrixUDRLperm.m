% The Green's matrix, which we derive from:
%Z = SUM(det(Md)det(Mup)), Gup and Gdown are Md and Mup

% expm(-deltaTau*Kmu) matrix is precomputed...
function [Gup,Gdown] = GreenMatrixUDRLperm(A,S, lambda,Lperm)
    
    [N,L] = size(S);
    
    v_1up = createV(S,Lperm(1),lambda);
    B1 = A*expm(v_1up);
    [Ql,R1, P1] = qr(B1); 
    Dl = diag(diag(R1));
    Tl = Dl^-1*R1*P1.';
    Tproduct = Tl;
    for l = 2:L
       lindex = Lperm(l);
       v_lup = createV(S,lindex,lambda);
       B_l = A*expm(v_lup);
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
    v_ldown = -createV(S,Lperm(1),lambda);
    B1d = A*expm(v_ldown);

    [Ql, R1, P1] = qr(B1d); 
    Dl = diag(diag(R1));
    Tl = Dl^-1*R1*P1.';
    Tproduct = Tl;
    for l = 2:L
        
       lindex = Lperm(l);
       v_ldown = -1*createV(S,lindex,lambda);
       B_l = A*expm(v_ldown);
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