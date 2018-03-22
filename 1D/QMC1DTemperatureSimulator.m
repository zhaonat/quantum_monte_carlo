%%QMC 1D simulator
addpath(genpath('D:\Nathan\Documents\StanfordYearOne\DevereauxGroup\MatlabQMC'))
MSGID = 'MATLAB:nearlySingularMatrix'
warning('off', MSGID)

%% ================ Simulation Parameters =====================
N = 6;
mu = 0;
t = 0.5;
iter = 1500;
burnin = 500;
%% ==========================================================

figure()
Ustrengths = [0.5,1,2,4,6,8,10,12];
Ustrengths = [0.5, 8];
%low temperature limit: roughly T = 1 to 10;
% Generate Temp: (numpoints, L, Umax, t); need U and t to enforce
% convergence
[imaginaryTimeSteps,Temps, TimeSlices] = GenerateTemperatureRange(10,20, max(Ustrengths),t,-0.75);
cmap = hsv(10)
counter = 1;
for U = Ustrengths
    MeanMoments = []; Errors = [];
    L = TimeSlices(counter);
    for deltaTau = imaginaryTimeSteps
        deltaTau
        %do several hundred iterations without measurements
        % do several thousand with measurements
        [Gup, Gdown, magneticMoms,S] = runQMC1D(N,U,mu,t,L,deltaTau,iter);
        %magneticMoms is an array containing all measurements for the
        %moment for each site
        
        %calculate ensemble properties
        SiteMeasurements = runQMCMeasurements(magneticMoms);
        SiteMoments = SiteMeasurements(:,1); SiteErrors = SiteMeasurements(:,2);
        %now condense site properties to single value using a mean
        MeanMoments = [MeanMoments, mean(SiteMoments)];
        Errors = [Errors, mean(SiteErrors)];

    end
    semilogx(Temps, MeanMoments, '.-', 'markersize', 20, 'color', rand(1,3))
    hold on;
    errorbar(Temps, MeanMoments, Errors);
    counter = counter+1;
end
%% reaching higher temperatures, can spread out Temps probed;

ylim([0.45,1.0])

grid()