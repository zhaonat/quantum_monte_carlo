%% Basic Sherman Morrison Formula

A = rand(10,10);

u = rand(10,1);
v = rand(10,1);

%%(A+uv^t)^-1 = A^-1 + A^-1*uv^T*A^-1/(1 + v^TA^-1*u)
LHS = (A + u*v.')^-1;
RHS = A^-1 - (A^-1*(u*v.')*A^-1)/(1 + v.'*(A^-1*u));

norm(LHS-RHS)
