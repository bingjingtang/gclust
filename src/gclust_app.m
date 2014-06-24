function [Xorigin, What, Hhat] = gclust_app(Xinput,innDim,nmfmethod)

    Xorigin = Xinput;
    Mhat = Xinput/diag(sum(Xinput));

    if nmfmethod==1 
        try
            Khat = sort(FastSepNMF(Mhat,innDim));
            What = Mhat(:,Khat);
            Hhat = nnlsHALSupdt(Mhat,Mhat(:,Khat));
        catch
            disp('FastSepNMF problem')
        end
    elseif nmfmethod==2
        try
            Khat = sort(FastConicalHull(Mhat,innDim));
            What = Mhat(:,Khat);
            Hhat = nnlsHALSupdt(Mhat,Mhat(:,Khat));
        catch
            disp('FastConicalHull Problem')
        end
    else 
        try
            [What,Hhat] = nnmf(Mhat,innDim,'algorithm','als','option',statset('TolX',1e-16,'TolFun',1e-16,'MaxIter',1e12));
        catch
            disp('nnmf problem')
        end
    end
    Hhat = diag(sum(What)) * Hhat ;
    Hhat = Hhat/diag(sum(Hhat));
    What = What/diag(sum(What));
end