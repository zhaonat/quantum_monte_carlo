%%QMC 1D simulator
close all

addpath('D:\Nathan\Documents\StanfordYearOne\DevereauxGroup\MatlabQMC')
N = 10;
mu = 0;
L = 10
deltaTau = 0.1
t = 1;
U = 2;

figure()
MeanMoments = []; Errors = [];
Nu = []; Nd = [];
InteractionStrengths = 0.01:0.2:3
path
[imaginaryTimeSteps,Temps, TimeSlices] = GenerateTemperatureRange(50,20, U,t,-0.75);

for deltaTau = imaginaryTimeSteps
    deltaTau
    iter = 1000;
    %do several hundred iterations without measurements
    % do several thousand with measurements
    [Gup, Gdown, magneticMoms,S] = runQMC1D(N,U,mu,t,L,deltaTau,iter);
    imagesc(Gup) %% it true that Gup = -Gdown?

    
    %calculate ensemble properties
    SiteMeasurements = runQMCMeasurements(magneticMoms);
    SiteMoments = SiteMeasurements(:,1); SiteErrors = SiteMeasurements(:,2);
    %now condense site properties to single value using a mean
    MeanMoments = [MeanMoments, mean(SiteMoments)];
    Errors = [Errors, mean(SiteErrors)]; 
end

semilogx(Temps, MeanMoments, '.-', 'markersize', 20, 'color', rand(1,3))
xlabel('temperatures')
ylabel('<mz>^2/');
grid()
hold on;
errorbar(Temps, MeanMoments, Errors);
% plot(Temps, Nu, '.', 'markersize', 10)
% plot(Temps, Nd, '.', 'markersize', 10)