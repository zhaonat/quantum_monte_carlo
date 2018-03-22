%%QMC 2D simulator
%turn off warnings for the poorly conditioned case
MSGID = 'MATLAB:nearlySingularMatrix'
warning('off', MSGID)

Nx = 6; Ny = Nx;
mu = 0; %half-filling
t = 1;
figure()

counter = 1;
iter = 500;
Ustrengths = [0.5, 1, 2];
% Generate Temp: (numpoints, L, Umax, t); need U and t to enforce
% convergence
errThresh = 0.1;
[imaginaryTimeSteps,Temps, TimeSlices] = ...
    GenerateTemperatureRange(15,10, max(Ustrengths),t,-1.1, errThresh) ;
legendInfo = cell(2,1);

for U = Ustrengths
    MeanMoments = [];
    for deltaTau = imaginaryTimeSteps
        L = TimeSlices(counter);
        
        cluster = 1;
        if(deltaTau < 0.05); cluster = 2; end;
        [Gup, Gdown, magneticMoms] = ...
            runQMC2D(Nx, Ny, U, mu, t, L, deltaTau, iter, cluster); 
        MeanMoments = [MeanMoments mean(magneticMoms)]

    end
    semilogx(Temps, MeanMoments, '.-', 'markersize', 20, 'color', rand(1,3))
    hold on;
    legendInfo{counter} = ['X = ' num2str(U)];
    counter = counter+1;
end

legend(legendInfo);
ylim([0.45,1.0])
grid();
xlabel('Temperatures (K)')
ylabel('Value of Moment')

