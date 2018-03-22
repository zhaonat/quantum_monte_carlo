function AutoCorr = SiteCorrelations(Gup, Gdown, Nx, Ny)
    N = length(Gup);
    AutoCorr = zeros(Nx, Ny);
    for i = 1:N
       for j = 1:N
           AutoCorr(i,j) = -Gup(j,i)*Gdown(i,j);
       end
    end
    
end