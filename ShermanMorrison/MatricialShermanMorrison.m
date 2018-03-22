function inverseAns = MatricialShermanMorrison(A, u, v)
    Ainv = A^-1;
    inverseAns = Ainv - (Ainv*(u*v.')*Ainv)/(1 + v.'*(Ainv*u));
end