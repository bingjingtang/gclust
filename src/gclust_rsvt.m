function [Xorigin, What, Hhat] = gclust_rsvt(Xinput,innDim, maxsvt)

    Xorigin = Xinput;
    Xinput = Xinput/diag(sum(Xinput));

    [U,S,V] = svds(Xinput,innDim);
    
    for itr = 1:maxsvt
        X = U*S*V';
        X(X<0) = 0;
        [UU,SS,VV] = svds(X,innDim);
        if norm(UU-U)+norm(VV-V)+norm(SS-S) < 1e-12
            break
        end
        U=UU;
        V=VV;
        S=SS;
    end

    XX = X/diag(sum(X));

    try
        [WW,HH] = nnmf(XX,innDim,'option',statset('TolX',1e-16,'TolFun',1e-16,'MaxIter',1e12));
    catch
        disp('nnmf problem')
    end
    
    Hhat = diag(sum(WW)) * HH ;
    Hhat = Hhat/diag(sum(Hhat));
    What = WW/diag(sum(WW));
end

