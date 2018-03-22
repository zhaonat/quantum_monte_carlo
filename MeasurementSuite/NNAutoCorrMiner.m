function NNAuto = NNAutoCorrMiner(AutoCorr)
    NNAuto = [];
    for i = 1:length(AutoCorr)
        NNR = i+1; NNL = i-1;
        if(i == length(AutoCorr)); NNR = 1; end;
        if(i == 1); NNL = length(AutoCorr);end;
        NNAuto = [NNAuto; AutoCorr(i,NNR), AutoCorr(i,NNL)];
        
    end
end