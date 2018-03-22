%% Matrix Decomposition Stabilization Analysis

addpath(genpath('D:\Nathan\Documents\StanfordYearOne\DevereauxGroup\MatlabQMC'))
clear;
cluster = 1;
N = 20;
% Ugrid = [0.5, 1, 2, 3, 5, 8, 10, 12]
U = 10;  %U = 0 means that all site changes are immediately accepted...interaction does not penalize
     % hopping?
%check the tutorial for half-filling
mu = U/2;%a nonzero mu fills the entire diagonal of Kmu matrix
t = 1; %t fills all the immediate off diagonals of Kmu

L = 10;
deltaTau = 0.2;
Beta = L*deltaTau; T = 1/Beta;
lambda = (acosh(exp(U*deltaTau/2))); 

S = createS(N,L);
%imag1esc(S)
Kmu = KEMatrix(N, t);  %mu*eye(N);
ExpmK = expm(-deltaTau*Kmu);

%% guts of the fast qr code
[N,L] = size(S);

v_lup1 = createV(S,L,lambda);
v_lup2 = createV(S,L-1,lambda);
v_lup3 = createV(S,L-2,lambda);

% B1 = ExpmK*expm(v_lup3)*ExpmK*expm(v_lup2)*ExpmK*expm(v_lup1);
B1 = eye(N); Lperm = L:-1:1;
cluster = 2;
for i = 1:cluster
    v_lup = createV(S, Lperm(i), lambda);
    B1 = ExpmK*expm(v_lup)*B1;
end

[Ql,R1, P1] = qr(B1); 
Dl = diag(diag(R1));
Tl = Dl^-1*R1*P1.';
Tproduct = Tl;

for l = cluster+1:cluster:L
   lindex = Lperm(l);
   B_l = eye(N);
   for i = l:l+cluster-1
        X = [num2str(l),', ',num2str(i)];
        display(X)
        v_lup = createV(S, Lperm(i), lambda);
        B_l = ExpmK*expm(v_lup)*B_l; %%order of multiplication matters
   end
  
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

[Gup, Gdown] = GreenMatrixLPerm(ExpmK, S, lambda,Lperm);
[Gup2, Gdown2] = GreenMatrix(deltaTau, Kmu, S, lambda);

disp('Clustered vs Full Green Matrix With L Perm')
norm(G-Gup)
norm(G-Gup2)

%use GreenMatrixFunc
Lperm = L:-1:1;
tic
[Gudru, Gudrd] = GreenMatrixUDRLperm(ExpmK,S, lambda, Lperm);
toc
norm(Gudru-Gup)
%use clustered Green MatrixFunc
tic
cluster = 10;
[Gudruf, Gudrdf] = GreenMatrixUDRFastLperm(ExpmK,S, lambda, Lperm, cluster);
toc


norm(Gudruf-Gup)
disp('clustered L-perm vs GreenMatrix')
norm(Gudruf-Gup2)
disp('clustered L-perm vs GreenMatrixLPerm')
norm(Gudrdf-Gdown)