
function [Gup, Gdown, magneticMom] = runQMC2D(Nx, Ny, U, mu, t, L, deltaTau, iter, cluster)
   %% Parameters
    if nargin < 9
        cluster =  1;
    end
    Temperature = 1/(deltaTau*L);
    UDRRate = 1;

    
    lambda = (acosh(exp(abs(U)*deltaTau/2)));
    %Creating the HubbardStratonovich Field
    S = HS2DField(Nx, Ny, L);
    N= Nx*Ny;
    Kmu = KE2DMatrix(Nx, Ny, t) + mu*eye(N);
    A = expm(-deltaTau*Kmu);
    
    Lperm = 1:L;
    %Create the Green's matrix
    [Gup, Gdown] = GreenMatrixUDRFastLperm(A, S, lambda,Lperm, cluster);
    thresh = round(iter/2);
    
    magneticMom = 0; nup =0; ndown = 0; ndouble = 0;
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
              if(r <d)
                  S0 = S;
                  S(i,l) = -S(i,l);
                  % update Green's matrices with new S
                  [Gup, Gdown] = ShermanMorrisonUpdate(Gup, Gdown, S, i, l, S0, lambda);
                  
              end

           end
            
           %% WRAP/RECOMPUTE Green's Functions
           if(mod(iteration, 1) == 0 && mod(l,UDRRate) == 0)
               [Gup, Gdown] = GreenMatrixUDRFastLperm(A,S,lambda,Lperm, cluster);
           else
               v_lup = createV(S,l,lambda); %update V for wrapping
               v_ldown = -v_lup;
               Gup = (A*expm(v_lup))*Gup* (A*expm(v_lup))^-1;
               Gdown = (A*expm(v_ldown))*Gdown* (A*expm(v_ldown))^-1;
           end
           
        end
        %% run measurements after certain number of iterations
        if(iteration > thresh)
            m1 = 1 - (diag(Gup)); m2 = 1-(diag(Gdown)); 
            m3 = (1 - (diag(Gup))).*(1-(diag(Gdown)))
            nup = nup +1 - (diag(Gup)); 
            ndown = ndown + 1-(diag(Gdown));
            ndouble = ndouble + (1 - (diag(Gup))).*(1-(diag(Gdown)))
            magneticMom = magneticMom + m1 + m2 - 2*m3;
        end

    end %% end of iterations
    %%normalize measurements
    normalize = iter-thresh;
    magneticMom = magneticMom/normalize;
    nup = nup/normalize;
    ndown = ndown/normalize;
    ndouble = ndouble/normalize;
    
end
