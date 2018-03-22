%%QMC 2D simulator
%turn off warnings for the poorly conditioned case
MSGID = 'MATLAB:nearlySingularMatrix'
warning('off', MSGID)

Nx = 6; Ny = Nx;
mu = 0; %half-filling
t = 1;
figure()
iter = 500;

%%================SELECT U PARAMETER =====================%
U = 0.5
%==========================================================%

% Generate Temp: (numpoints, L, Umax, t); need U and t to enforce
% convergence
[deltaTaus,Temps, TimeSlices] = ...
    GenerateTemperatureRange(20,10, U,t,-1.4,0.1) ;legendInfo = cell(2,1);
MeanMoments = [];
counter = 1;
for deltaTau = deltaTaus
    L = TimeSlices(counter);
    cluster = 1;
    if(deltaTau < 0.05); cluster = 2; end;
    [Gup, Gdown, magneticMoms] = ...
        runQMC2D(Nx, Ny, U, mu, t, L, deltaTau, iter,cluster); 
    MeanMoments = [MeanMoments mean(magneticMoms)]
    counter = counter+1;

end

semilogx(Temps, MeanMoments, '.-', 'markersize', 20, 'color', rand(1,3))
legend(strcat('U=',num2str(U)));
ylim([0.45,1.0])
grid();
xlabel('Temperatures (K)')
ylabel('Value of Moment')

