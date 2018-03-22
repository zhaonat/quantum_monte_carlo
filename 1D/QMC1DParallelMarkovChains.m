%%QMC 1D simulator

addpath('D:\Nathan\Documents\StanfordYearOne\DevereauxGroup\MatlabQMC')
N = 10;
mu = 0;
L = 10
deltaTau = 0.1
t = 1;
U = 2;
numChains = 8;

    
%% Variable initialization for each Markov Chain
MeanMoments = []; Errors = [];
Nu = []; Nd = [];
path
[imaginaryTimeSteps,Temps, TimeSlices] = ...
    GenerateTemperatureRange(20,20, U,t,-0.75, 1, 0.01);
    
 for deltaTau = imaginaryTimeSteps
    MarkovChains = [cell(1),cell(1),cell(1),cell(1)];
    parfor i = 1:numChains
        deltaTau
        iter = 200;
        %do several hundred iterations without measurements
        % do several thousand with measurements
        [Gup, Gdown, magneticMoms,S] = runQMC1D(N,U,mu,t,L,deltaTau,iter);
        imagesc(Gup) %% it true that Gup = -Gdown?

        %calculate ensemble properties
        bin = 20;
        [SiteMeasurements, ClustMeasurements] = runQMCMeasurements(magneticMoms, bin);
        SiteMoments = ClustMeasurements(:,1); SiteErrors = ClustMeasurements(:,2);
        %now condense site properties to single value using a mean

        MarkovChains{i} = SiteMeasurements;

    end
    
    %% aggregate the parallel markov chain measures
    siteMeasure = 0;
    for j = 1:numChains
        siteMeasure =siteMeasure + MarkovChains{j};
    end
    MeanMoments = [MeanMoments, mean(siteMeasure(:,1))/numChains];
    Errors = [Errors, std(siteMeasure(:,2))/numChains]; %divided by sqrt(n) in the function 
    
 end

%% PLOT FIGURES
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

%% 