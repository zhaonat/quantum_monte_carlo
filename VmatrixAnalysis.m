%% analysis of the V matrix (decoupled interaction term in Hubbard)
N = 20;
L = 5; lambda = 0.1
S = createS(N, L)

for i = 1:L
   figure()
   Vup = createV(S, i, 1, lambda);
   Vdown = -Vup;
   subplot(2,1,1)
   imagesc(Vdown)
   subplot(2,1,2)
   imagesc(Vup)
end