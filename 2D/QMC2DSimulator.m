%%QMC 2D simulator
%turn off warnings for the poorly conditioned case
MSGID = 'MATLAB:nearlySingularMatrix'
warning('off', MSGID)
dir = 'D:\Nathan\Documents\StanfordYearOne\DevereauxGroup\MatlabQMC\'
addpath(genpath(dir));
Nx = 6; Ny = Nx;
mu = 0; %half-filling
t = 1;
figure()
counter = 1;
iter = 1000;
Ustrengths = [0.5,1,2,4,6,8,10,12];

% Generate Temp: (numpoints, L, Umax, t); need U and t to enforce
% convergence
[imaginaryTimeSteps,Temps, TimeSlices] = ...
    GenerateTemperatureRange(1000,10, max(Ustrengths),t,-2) ;
figure()
semilogx(Temps, TimeSlices, 'linewidth', 2)
xlabel('Temps')
ylabel('L (# of Time Slices)')
title('Scaling of Time Slices To Satisfy Trotter Error Reqs')
legendInfo = cell(2,1);
MomentsMatrix = [];
for U = Ustrengths
    MeanMoments = [];
    for deltaTau = imaginaryTimeSteps
        L = TimeSlices(counter);
        
        cluster = 2;
        [Gup, Gdown, magneticMoms] = ...
            runQMC2D(Nx, Ny, U, mu, t, L, deltaTau, iter,cluster); 
        MeanMoments = [MeanMoments mean(magneticMoms)]

    end
    semilogx(Temps, MeanMoments, '.-', 'markersize', 20, 'color', rand(1,3))
    hold on;
    legendInfo{counter} = ['X = ' num2str(U)];
    counter = counter+1;
    MomentsMatrix = [MomentsMatrix; MeanMoments];
end
%csvwrite(strcat(dir,'\MeanMoments.csv'), MomentsMatrix);

legend(legendInfo);
ylim([0.45,1.0])
grid();
xlabel('Temperatures (K)')
ylabel('Value of Moment')

