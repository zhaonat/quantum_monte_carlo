% The Green's matrix, which we derive from:
%Z = SUM(det(Md)det(Mup)), Gup and Gdown are Md and Mup
% FAST UDR
% expm(-deltaTau*Kmu) matrix is precomputed...
function [Gup,Gdown] = GreenMatrixUDRFastLperm(ExpmK,S, lambda,Lperm, cluster)
    
    [N,L] = size(S);
    B1 = eye(N); 
    for i = 1:cluster
        v_lup = createV(S, Lperm(i), lambda);
        B1 = ExpmK*expm(v_lup)*B1;
    end
    [Ql,R1, P1] = qr(B1); 
    Dl = diag(diag(R1));
    Tl = Dl^-1*R1*P1.';
    for l = cluster+1:cluster:L

        B_l = eye(N);
        for i = l:l+cluster-1
            if(i>L); continue; end
%             X = [num2str(l),', ',num2str(i)];
%             display(X)
            v_lup = createV(S, Lperm(i), lambda);
            B_l = ExpmK*expm(v_lup)*B_l; %%order of multiplication matters
        end
        
       C_l = (B_l*Ql)*Dl; 
       [Ql,R,P] = qr(C_l); %q*r*p.' = C...the original matrix, got to transpose P
       Dl = diag(diag(R));
       Tl = Dl^-1 * R * P.'*Tl;
    end   
    %we have to get the inverse of the Green's matrix
    Db = DecompDb(Dl);
    Ds = DecompDs(Dl);
    % Db*Ds = Dl
    %norm(Db*Ds-Dl)
    Gup = (Db^-1*Ql.' + Ds*Tl)^-1*(Db^-1*Ql.');

    %% Now do the same for Gdown
    
    B1 = eye(N); 
    for i = 1:cluster
        v_ldown = -createV(S, Lperm(i), lambda);
        B1 = ExpmK*expm(v_ldown)*B1;
    end
    [Ql,R1, P1] = qr(B1); 
    Dl = diag(diag(R1));
    Tl = Dl^-1*R1*P1.';
    for l = cluster+1:cluster:L

        B_l = eye(N);
        for i = l:l+cluster-1
            if(i>L); continue; end
%             X = [num2str(l),', ',num2str(i)];
%             display(X)
            v_ldown = -createV(S, Lperm(i), lambda);
            B_l = ExpmK*expm(v_ldown)*B_l; %%order of multiplication matters
        end
        
       C_l = (B_l*Ql)*Dl; 
       [Ql,R,P] = qr(C_l); %q*r*p.' = C...the original matrix, got to transpose P
       Dl = diag(diag(R));
       Tl = Dl^-1 * R * P.'*Tl;
    end   
    %we have to get the inverse of the Green's matrix
    Db = DecompDb(Dl);
    Ds = DecompDs(Dl);
    Gdown = (Db^-1*Ql.' + Ds*Tl)^-1*(Db^-1*Ql.');
    
end