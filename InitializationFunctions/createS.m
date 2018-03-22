%% generate the hubbard stratonovitch field

function S = createS(N,L)

%     S = randi([0,1], N,L);
    S=random('Uniform',-1,1,[N,L]);

    %convert all 0's to -1
%     for idx = 1:numel(S)
%         if S(idx) == 0
%            S(idx) = -1 ; %%why are elements plus/minus 1
%         end
%     end
    
end