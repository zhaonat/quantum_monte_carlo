function S2D = HS2DField(Nx, Ny, L)
    %each slice is Nx*Ny, indexing positions in natural ordering
    %then we have L slices in imaginary time
    N = Nx*Ny;
    S = randi([0,1], N,L);
    %convert all 0's to -1
    for idx = 1:numel(S)
        if S(idx) == 0
           S(idx) = -1 ;
        end
    end
    S2D = S;
    
end