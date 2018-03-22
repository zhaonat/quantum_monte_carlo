%% create a matrix that has dimension Nx*Ny;
% this matrix should theoretically be diagonal, as the onsite interaction
% only affects the site itself, not reaching out to others

function FC = ImpurityInteractionMatrix(Nx, Ny, V_k)
    FC = V_k*eye(Nx*Ny);
end