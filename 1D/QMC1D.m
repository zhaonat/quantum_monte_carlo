close all
[a, MSGID] = lastwarn();
warning('off', MSGID)
addpath(genpath('D:\Nathan\Documents\StanfordYearOne\DevereauxGroup\MatlabQMC'))
N = 20;
U = 10;  %U = 0 means that all site changes are immediately accepted...interaction does not penalize
     % hopping?
%check the tutorial for half-filling
mu = 0;
t = 1; %t fills all the immediate off diagonals of Kmu
%beta = deltaTau*L

% Set L first; L > 2, then deltaTau, then Temp;  
L = 20;
deltaTau = 0.1;
Beta = L*deltaTau; T = 1/Beta;

lambda = (acosh(exp(abs(U)*deltaTau/2))); 

S = createS(N,L);

%imagesc(S)
Kmu = KEMatrix(N, t);  %mu*eye(N);
A = expm(-deltaTau*Kmu);
Lperm = L:-1:1;
[Gup, Gdown] = GreenMatrixUDRLperm(A,S, lambda,Lperm);
Gu0 = Gup; Gud = Gdown;

iter = 50; thresh = iter/2; %threshold for measuring shit
v_lup = 0; v_ldown = 0;
magneticMom = []; nup = []; ndown = []; ndouble = [];
D = ones(N*L*iter,3); % analyze the probabilistic weight ratios
numAccepts = 0; counter = 0;

%%note that we still get numerical error accumulation, so we should
for iteration = 1:iter
    iteration
    
    for l = 1:L %based on how we created the green's matrix, the Lth index is at the end!!!
       %should we initialize the green's matrix elements for slice i? 
       for i = 1:N
          stationary = S(i,l);
          % G ~ <cicj*>           
          dup = 1 + (1 - Gup(i,i))*(exp(-2*lambda*stationary)-1); 
          ddown = 1 + (1 - Gdown(i,i))*(exp(2*lambda*stationary)-1);
          d = (dup*ddown) ;
          D(counter+1, :) = [dup, ddown, d];
         
          r = rand(1);
          counter = counter+1;
          if(r < min(1,d)) %remember, we want to favor higher weights
              S0 = S;
              S(i,l) = -S(i,l);
              % update Green's matrices, Sherman morrison requires that the
              % lth term be on the VERY RIGHT
              [Gup, Gdown] = ShermanMorrisonII(Gup, Gdown, S, i, l, S0, lambda);
              numAccepts  = numAccepts + 1;

          else
              continue; %no updates to the Green's matrix
          end
          
       end %end of N loop through sites

       %% wrapping the green's function move next exp(k)exp(v_l) to front
       % or recompute the Green's function with correct L ordering
       
       Lperm = circshift(Lperm,1);
       if(mod(iteration, 1) == 0 && mod(l,5) == 0) %% inside L loop;
           [Gup, Gdown] = GreenMatrixUDRLperm(A,S,lambda,Lperm);
       else
           v_lup = createV(S,l,lambda); %update V for wrapping
           v_ldown = -v_lup;
           %this update introduces a lot of error
           Blup = A*expm(v_lup);
           Bldown = A*expm(v_ldown);
           Gup = (Blup*Gup)*Blup^-1;
           Gdown = (Bldown*Gdown)*Bldown^-1;
       end
    end %end o fL loop

    %% Calculate Ensemble Quantities
    if( iteration > thresh ) %start doing measurements
        m1 = 1 - (diag(Gup)); m2 = 1-(diag(Gdown)); 
        m3 = (1 - (diag(Gup))).*(1-(diag(Gdown)));
        nup = [nup,m1]; 
        ndown = [ndown,m2];
        ndouble = [ndouble,m3];
        magneticMom = [magneticMom, m1 + m2 - 2*m3];
    end
    
end

figure()
plot(mean(magneticMom.'), 'linewidth', 2)
hold on;
plot(mean(magneticMom.')+std(magneticMom.')/sqrt(iter));
plot(mean(magneticMom.')-std(magneticMom.')/sqrt(iter));

figure()
plot(mean(nup.'))
hold on;
plot(mean(ndown.'))
plot(mean(ndouble.'))
legend('nup', 'down', 'ndown+nup')


