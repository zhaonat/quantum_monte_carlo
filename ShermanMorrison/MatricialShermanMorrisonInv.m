function inverseAns = MatricialShermanMorrisonInv(Ainv, u, v)
    inverseAns = Ainv - (Ainv*(u*v.')*Ainv)/(1 + v.'*(Ainv*u));
end