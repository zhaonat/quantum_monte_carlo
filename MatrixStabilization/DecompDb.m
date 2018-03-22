
%helper function for matrix decomposition stabilization
function Db = DecompDb(D_L)
    Db = zeros(length(D_L));
    for i = 1:length(D_L)
        d = diag(D_L);
        if(abs(d(i)) >1)
           Db(i,i) = (d(i));
        else
           Db(i,i) = 1;
        end
    end
end