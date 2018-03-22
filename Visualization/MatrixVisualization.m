%% Visualize Matrix
close all
A = KE2DMatrix(6,6,1);
imagesc(A);
title('Hopping Matrix for Conduction Electrons in 2d')

%% randomly initialized Hubbard Stratonovitch field
B = createS(36,20);
imagesc(B)
title('HS Field of Ising Spins')
xlabel('Imaginary Time Axis')
ylabel('Real-Space Lattice')
colormap('winter')
