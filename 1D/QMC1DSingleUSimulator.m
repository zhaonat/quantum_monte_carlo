%%QMC 1D simulator

addpath('D:\Nathan\Documents\StanfordYearOne\DevereauxGroup\MatlabQMC')
N = 10;
mu = 0;
L = 10
deltaTau = 0.1
t = 1;
U = 10;

figure()
MeanMoments = []; Errors = [];
Nu = []; Nd = [];
path
[imaginaryTimeSteps,Temps, TimeSlices] = GenerateTemperatureRange(20,20, U,t,-0.75, 1, 0.01);

for deltaTau = imaginaryTimeSteps
    deltaTau
    iter = 300;
    %do several hundred iterations without measurements
    % do several thousand with measurements
    [Gup, Gdown, magneticMoms,S] = runQMC1D(N,U,mu,t,L,deltaTau,iter);
    imagesc(Gup) %% it true that Gup = -Gdown?

    
    %calculate ensemble properties
    bin = 20;
    [SiteMeasurements, ClustMeasurements] = runQMCMeasurements(magneticMoms, bin);
    SiteMoments = SiteMeasurements(:,1); SiteErrors = ClustMeasurements(:,2);
    %now condense site properties to single value using a mean
    MeanMoments = [MeanMoments, mean(SiteMoments)];
    Errors = [Errors, mean(SiteErrors)]; %divided by sqrt(n) in the function 
    
end
fig = figure()
semilogx(Temps, MeanMoments, '.-', 'markersize', 20, 'color', rand(1,3))
xlabel('temperatures')
ylabel('<mz>^2/');
grid()
hold on;
errorbar(Temps, MeanMoments, Errors);
% plot(Temps, Nu, '.', 'markersize', 10)
% plot(Temps, Nd, '.', 'markersize', 10)
saveas(fig, strcat('Single U=',num2str(U)))