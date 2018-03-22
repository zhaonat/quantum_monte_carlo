%% 2D Quantum Monte Carlo Determinant
close all
clear all
addpath(genpath('D:\Nathan\Documents\StanfordYearOne\DevereauxGroup\MatlabQMC'))


%% Parameters
Nx = 6; Ny = 6;

%% Ugrid = [0.5, 1, 2, 3, 5, 8, 10, 12]
U = 4
mu = 0;
t = 1; L = 10;
deltaTau = 0.1;
lambda = acosh(exp(abs(U)*deltaTau/2));
beta = deltaTau*L;
T = 1/beta;

%Creating the HubbardStratonovich Field
% index as (1,1),(2,1),(3,1)...(1,2),(2,2)
S = HS2DField(Nx, Ny, L);

%it still seems like we can use the original createV function to generate
%the V's for L different time slices (the V's will remain diagonal)
N = Nx*Ny;
Kmu = KE2DMatrix(Nx, Ny, t) + mu*eye(N);
Lperm = 1:L;
A = expm(-deltaTau*Kmu);
%Create the Green's matrix
[Gup, Gdown] = GreenMatrixLPerm(A, S, lambda, Lperm);
v_lup = 0; v_ldown = 0;
iter = 100; thresh = round(iter/2);
UDRRate = 2;
nup = 0; ndown = 0; ndouble = 0; magneticMom = 0;
Autocorr = zeros(N, N);
for iteration = 1:iter
    iteration
    for l = L:-1:1
       Lperm = circshift(Lperm,1);

       for i = 1:N
          stationary = S(i,l);
          dup = 1 + (1 - Gup(i,i))*(exp(-2*lambda*stationary)-1);
          ddown = 1 + (1 - Gdown(i,i))*(exp(+2*lambda*stationary)-1);
         
          d = dup*ddown;
          r = rand();
          if(r < d)
              S0 = S;
              S(i,l) = -S(i,l);
              % update Green's matrices with new S
              [Gup, Gdown] = ShermanMorrisonUpdate(Gup,Gdown, S, i, l,S0, lambda);

          end
          
       end
       
       %% WRAP/RECOMPUTE Green's Functions
       if(mod(iteration, 1) == 0 && mod(l,1) == 0)
           [Gup, Gdown] = GreenMatrixUDRFastLperm(A,S,lambda,Lperm,5);
       else
           v_lup = createV(S,l,lambda); %update V for wrapping
           v_ldown = -v_lup;
           Gup2 = (A*expm(v_lup))*Gup* (A*expm(v_lup))^-1;
           Gdown2 = (A*expm(v_ldown))*Gdown* (A*expm(v_ldown))^-1;
       end

    end
    %% Measure Error
    
    %% Calculate Ensemble Quantities
    if( iteration > thresh ) %start doing measurements
        m1 = 1 - (diag(Gup)); m2 = 1-(diag(Gdown)); 
        m3 = (1 - (diag(Gup))).*(1-(diag(Gdown)));
        nup = nup +1 - (diag(Gup)); 
        ndown = ndown + 1-(diag(Gdown));
        ndouble = ndouble + (1 - (diag(Gup))).*(1-(diag(Gdown)));
        magneticMom = magneticMom + m1 + m2 - 2*m3;
        Autocorr = Autocorr - Gup.'*Gdown;

    end
  
end
normalization = iter-thresh
figure()
plot(magneticMom/normalization);
figure();
plot(nup/normalization)
hold on
plot(ndown/normalization)
plot(ndouble/normalization)

