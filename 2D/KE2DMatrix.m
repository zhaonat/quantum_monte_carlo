
%%this function should not be used for large grids
function KE2D = KE2DMatrix(Nx, Ny, t)
   
    KE2D = [];
    for i = 1:Nx
        for j = 1:Ny
           lattice = zeros(Nx, Ny);
           lattice(j,i) = 1
           lattice = lattice+circshift(lattice, 1, 1) + circshift(lattice, -1,1)+...
               circshift(lattice, 1,2) + circshift(lattice, -1,2)
           lattice(j,i) = 0; %% no self interaction hopping or mu terms
           row = t*reshape(lattice, 1,prod(size(lattice)))
           KE2D = [KE2D; row]
        end

    end
    
end

