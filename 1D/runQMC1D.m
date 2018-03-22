
function [Gup, Gdown, magneticMoms, S] = runQMC1D(N, U, mu, t, L, deltaTau, iter, burnin)
    if nargin < 7
        iter =  10;
    end
    if(nargin < 8)
       burnin = ceil(iter/2); 
    end
    Temperature = 1/(deltaTau*L);
    UDRRate = 5;
    if(Temperature < 1)
        UDRRate = 1;
    end
    magneticMoms = []; nup = []; ndown = []; ndouble = [];
    Lperm = 1:L;
    %%this formula comes from the discrete HS transformation
    lambda = acosh(exp(abs(U)*deltaTau/2));
    
    S = createS(N,L);
    %Kmu holds all non-interacting terms so chemical potential mu and KE
    Kmu = KEMatrix(N, t) + mu*eye(N);
    A = expm(-deltaTau*Kmu);
    [Gup, Gdown] = GreenMatrixUDRLperm(A, S, lambda, Lperm);

    for iteration = 1:iter
        iteration
        %for every time slice in the S-aux field 
        for l = L:-1:1
           Lperm = circshift(Lperm,1);
           %for every physical domain site
           for i = 1:N
              
              %% why do we only target diagonal elements?
              dup = 1 + (1 - Gup(i,i))*(exp(-2*lambda*S(i,l))-1);
              ddown = 1 + (1 - Gdown(i,i))*(exp(+2*lambda*S(i,l))-1);
              
              d = dup*ddown;
              r = rand();
              if(r < min(1,d)) %sometimes d is greater than 1, in which case we automatically accept 
                  %change S at the given site only if the dup*ddown isn't
                  %within range
                  S0 = S;
                  S(i,l) = -S(i,l); % -> V depends on S explicitly
                  % update Green's matrices with new S
                  [Gup, Gdown] = ShermanMorrisonII(Gup, Gdown, S, i, l, S0, lambda);
                  
              end
            

           end %% end site loop;
           
           if(mod(iteration, 1) == 0 && mod(l,UDRRate) == 0)
               cluster = 5;
               [Gup, Gdown] = GreenMatrixUDRFastLperm(A,S,lambda,Lperm, cluster);
           else
               v_lup = createV(S,l,lambda); %update V for wrapping
               v_ldown = -v_lup;
               Gup = (A*expm(v_lup))*Gup* (A*expm(v_lup))^-1;
               Gdown = (A*expm(v_ldown))*Gdown* (A*expm(v_ldown))^-1;
           end
        end %% end L imaginary time loop
        %% Calculate Ensemble Quantities
        if(iteration > burnin)
            m1 = 1-abs(diag(Gup)); m2 = 1- abs(diag(Gdown));
            m3 = m1.*m2;
            nup = [nup; 1-abs(diag(Gup))]; 
            ndown = [ndown; 1-abs(diag(Gdown))];
            ndouble = [ndouble; (1-abs(diag(Gup))).*(1-abs(diag(Gdown)))];
            magneticMoms =  [magneticMoms, m1 + m2- 2*m3];
        end
  
    end %end sweeps of the HS fields

end
