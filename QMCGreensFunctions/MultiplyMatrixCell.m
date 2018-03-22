function A = MultiplyMatrixCell(B)

    product = B{1};
    for i = 2:length(B)
        product = product*B{i};
    end
    A = product;

end