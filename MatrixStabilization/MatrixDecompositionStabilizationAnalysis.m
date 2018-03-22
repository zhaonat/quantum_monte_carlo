%% Matrix Decomposition Stabilization Analysis

addpath(genpath('D:\Nathan\Documents\StanfordYearOne\DevereauxGroup\MatlabQMC'))
clear;

N = 20;
% Ugrid = [0.5, 1, 2, 3, 5, 8, 10, 12]
U =10;  %U = 0 means that all site changes are immediately accepted...interaction does not penalize
     % hopping?
%check the tutorial for half-filling
mu = U/2;%a nonzero mu fills the entire diagonal of Kmu matrix
t = 1; %t fills all the immediate off diagonals of Kmu

L = 10;
deltaTau = 0.001;
Beta = L*deltaTau; T = 1/Beta;
lambda = (acosh(exp(U*deltaTau/2))); 

S = createS(N,L);
%imag1esc(S)
Kmu = KEMatrix(N, t);  %mu*eye(N);
A = expm(-deltaTau*Kmu);
[N,L] = size(S);

U_up = eye(N);D_up = eye(N);
Aup = eye(N); Adown = eye(N);

v_1up = createV(S,L,lambda);
B1 = expm(-deltaTau*Kmu)*expm(v_1up);
Lperm = L:-1:1;
[Ql,R1, P1] = qr(B1); 
Dl = diag(diag(R1));
Tl = Dl^-1*R1*P1.';

for l = 2:L
   lindex = Lperm(l);
   v_lup = createV(S,lindex,lambda);
   B_l = expm(-deltaTau*Kmu)*expm(v_lup);
   C_l = (B_l*Ql)*Dl; 
   [Ql,R,P] = qr(C_l); %q*r*p.' = C...the original matrix, got to transpose P
   Dl = diag(diag(R));
   Tl = Dl^-1 * R * P.'*Tl;
   
end

Bl = Ql*R*Tl;
%we have to get the inverse of the Green's matrix
Matrix = eye(N) + Bl;
Db = DecompDb(Dl);
Ds = DecompDs(Dl);
% Db*Ds = Dl
%norm(Db*Ds-Dl)
G = (Db^-1*Ql.' + Ds*Tl)^-1*(Db^-1*Ql.');

%% comparison
Lperm = 1:L;

[Gup, Gdown] = GreenMatrixWithLPerm(A, S, lambda,Lperm);

[Gup2, Gdown2] = GreenMatrix(deltaTau, Kmu, S, lambda);

norm(G-Gup)
norm(G-Gup2)

%use GreenMatrixFunc
Lperm = L:-1:1;
tic
[Gudru, Gudrd] = GreenMatrixUDRLperm(A,S, lambda, Lperm);
toc
norm(Gudru-Gup)
%use clustered Green MatrixFunc
tic
cluster = 1;
[Gudruf, Gudrdf] = GreenMatrixUDRFastLperm(A,S, lambda, Lperm, cluster);
toc

norm(Gudruf-Gup)
norm(Gudrdf-Gdown)