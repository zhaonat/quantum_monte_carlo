%% ShermanMorrison Analyzer
%% Notes
% Should the formula break down for high trotter errors?
% my intuition is no, since the SM formula is supposed to be exact
close all
addpath(genpath('D:\Nathan\Documents\StanfordYearOne\DevereauxGroup\MatlabQMC'))

N = 10;
U = 10  

mu = U/2;%a nonzero mu fills the entire diagonal of Kmu matrix
t = 1 %t fills all the immediate off diagonals of Kmu
%beta = deltaTau*L
%T = 50; %high T -> magnetization is 0.5, but we observe that gup, gdown become exactly 0.5 
        %at high T...

% Set L first; L > 2, then deltaTau, then Temp;  
L = 10
deltaTau = 1
Beta = L*deltaTau; T = 1/Beta;

lambda = ((acosh(exp(abs(U)*deltaTau/2)))); 

S = randi([0,1], N,10);
%convert all 0's to -1
for idx = 1:numel(S)
    if S(idx) == 0
       S(idx) = -1 ;
    end
end
u = zeros(10,1);
v = zeros(10,1);

%imagesc(S)
Kmu = KEMatrix(N, t)  %mu*eye(N);

[Gup, Gdown] = GreenMatrix(deltaTau, Kmu,S, lambda);
Gu0 = Gup; Gud = Gdown;

%% change a value in HS field
i = 4; l = 10; %% no matter what index we flip the field value, the change in G
               % always occurs in the diagonal

S0 = S;
S(i,l) = -S(i,l);
v(i) = exp(-lambda*(S(i,l)-S0(i,l)))-1;
u(i) = 1;
% revaluate green's matrix: here we have to use the new S
[Gfcu, Gfcd] = GreenMatrix(deltaTau, Kmu, S, lambda);

%% test different sherman morrison functions
imagesc(Gfcu - Gup)
figure()
subplot(2,1,1)
imagesc(Gup)
[Gsmu, Gsmd, c_kup, b_jup] = ShermanMorrisonUpdate(Gup, Gdown, S, i, l, S0, lambda);
subplot(2,1,2)
imagesc(Gup)
[Gsm2u, Gsm2d] = ShermanMorrisonII(Gup, Gdown, S, i, l, S0, lambda);

%% use the matricial sherman morrison formula
Gmu = MatricialShermanMorrisonInv(Gup, u, v);
%Gmd = MatricialShermanMorrison(Gdown, u, v);

norm(Gfcu - Gmu)
%% error analysis in L2 norm
disp('gubernatis sherman morrison')
e1 = norm(Gfcu - Gsm2u)
e2 = norm(Gfcd - Gsm2d)

disp('original sherman morrison')
e3 = norm(Gfcu - Gsmu)
e4 = norm(Gfcd - Gsmd)

figure()
subplot(2,1,1)
imagesc(log(abs(Gfcu - Gsm2u)))
subplot(2,1,2)
imagesc(log(abs(Gfcu - Gsmu)))


