% data format should be a T*N matrix, where T is the number of iterations
% and N is the number of sites simulated

function [SiteStats,ClusterMeans] = runQMCMeasurements(data, bin)
    if(nargin <2)
        bin = 10; 
    end
    [T,N] = size(data)
    SiteStats = [];
    numClusts = N/bin;
    currentclust = [];
    ClusterMeans = [];
    for i = 1:N
       %get site averages
       %standard error = std/sqrt(T)
       SiteStats = [SiteStats; mean(data(:,i)), std(data(:,i))];
       currentclust = [currentclust;  mean(data(:,i)), std(data(:,i))];
       if(mod(i, numClusts) == 0)
          ClusterMeans = [ClusterMeans; mean(currentclust(:,1)),...
              std(currentclust(:,1)/sqrt(bin))];
          currentclust = [];
       end
    end
    
end