function Varray = V_spin_arr(S, sign, lambda)
    [N,L] = size(S);
    if(abs(sign) > 1)
        error('sign should be 1 or -1');
    end
    C = cell(1,L);
    Vmat = [];
    for i = 1:L
        V = sign*lambda*diag(S(:,i)); 
        C{i} = V;
    end
    Varray = C;
end
