%% create a map of where impurity sites are in PAM
%% this is a real spatial map, not a matrix
%% actually, what we really want is a map of indices (1,2,3,4...) to the spatial sites
%% corresponding to the impurities
function [impuritySites, impurityIndex] = MapImpurities(Nx, Ny)
    impuritySites = zeros(Nx, Ny);
    impurityIndex = [];
    counter = 1;
    for i = 1:Nx
        for j = 1:Ny
            if(mod(i,2) == 0 && mod(j,2) == 0)
               impuritySites(i,j) = 1;
               impurityIndex = [ImpurityIndex, counter]; 
               counter = counter+1;
            end
        end
    end
end