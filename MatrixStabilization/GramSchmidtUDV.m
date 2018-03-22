
%% Algorithm 25 from Gubernatis QMC methods
% V is upper triangular
% D is diagonal
% U is an orthonormal matrix; testing shows this is not true from output

function [U,D,V] = GramSchmidtUDV(U)
    N = length(U) % U has to be an NxN square matrix
    D = eye(N);;
    V = eye(N);
    for k = 1:N
        for i = 1:N
            D(k,k) = D(k,k) + U(i,k)*U(i,k);
        end
        D(k,k) = D(k,k)^0.5;
        for i =1:N
            U(i,k) = U(i,k)/D(k,k);
        end
        V(k,k) = 1;
        for j = k+1:N
            V(k,j) =0;
            for i = 1:N
                V(k,j) = V(k,j) + U(i,k)*U(i,j);
            end
            for i = 1:N
                U(i,j) = U(i,j) - V(k,j)*U(i,k);
            end
            V(k,j) = V(k,j)/D(k,k);
        end
    end
end