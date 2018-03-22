
%helper for matrix decomposition
function Ds = DecompDs(D_L)
    Ds = zeros(length(D_L));
    for i = 1:length(D_L)
        d = diag(D_L);
        if(abs(d(i)) <= 1)
           Ds(i,i) = d(i);
        else
           Ds(i,i) = 1;
        end
    end
end