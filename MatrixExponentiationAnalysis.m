%% Analysis of matrix exponentiation
N = 20;
mu = 0.1;
KE = KEMatrix(0.1, N, 1) + mu*eye(N);
[U, V] = eig(KE);

a = eye(N);
for i = 1:1000
   a = a +(1/factorial(i))*U*V^i*U^-1; 
end

figure()
imagesc(a);
figure()
imagesc(expm(KE)) 
%%verifying that expm(KE) is still a relatively sparse matrix
S = createS(N,N);
V = createV(S,1,1,0.1)
%does expm(V)expm(K) = expm(K)*expm(V)
