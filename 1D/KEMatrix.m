function KE = KEMatrix(N, t)
    I = eye(N);
    KE = circshift(I,-1)+circshift(I,1);
    KE = t*KE;
end