function B = MatrixMultCluster(S,ExpmK, cluster, Lperm, lambda)
    [N,L] = size(S)
    B = cell(ceil(L/cluster),1);
    for l = 1:cluster:L
       Bl = eye(N);
       for i = l:l+cluster-1
           lindex = Lperm(i);
           v = createV(S,lindex,lambda);
           Bl = Bl*ExpmK*expm(v);
       end
       index = ceil(l/cluster);
       B{index} = Bl;
    end
end