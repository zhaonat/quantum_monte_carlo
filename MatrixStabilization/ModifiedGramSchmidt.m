%simple gram-schmidt orthogonalization

function VectorMatrix = ModifiedGramSchmidt(VectorMatrix)
    M = length(VectorMatrix);
    for i = 1:M
        ai = VectorMatrix(:,i)/norm(VectorMatrix(:,i));
        VectorMatrix(:,i) = ai;
        for j = i+1:M
            aj = VectorMatrix(:,j)
            VectorMatrix(:,j) = VectorMatrix(:,j) - (ai.'*aj)*ai
        end
    end
    
end