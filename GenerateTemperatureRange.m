function [deltaTaus, Temps, TimeSlices] = ...
GenerateTemperatureRange(numPoints, Lstart,U,t, logLowLim, logHighLim, errThresh)
    if(nargin < 7)
         errThresh = 0.1
    end
    if(nargin < 6)
       LogHighLim = 2; 
    end
    deltaTaus = [];
    TimeSlices = [];
    %create a logarithmic scale
    Temps= logspace(logLowLim,logHighLim, numPoints);
    for T = Temps;
        L = Lstart;
        Beta = 1/T;

        deltaTau = Beta/L;
        
        %enforce condition that U*t*deltaTau^2 < 0.1;
        while(U*t*deltaTau^2 >= errThresh && deltaTau > 0.05)
            L = L+1;
            deltaTau = Beta/L;
        end
        TimeSlices = [TimeSlices L];
        deltaTaus = [deltaTaus deltaTau];
        
        
    end
end