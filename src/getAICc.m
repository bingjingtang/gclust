function [innDim, negloglik, penalty, AIC] = getAICc(Xorigin,What,Hhat)

    innDim = size(What,2);

    Phat = What * Hhat;
    Phat(Phat<1e-7) = 0;

    nparams = sum(What>0)/diag(sum(Xorigin * Hhat'));
    nparams = sum(nparams);
   
    zeroprob = find(Phat<1e-7);
    
    if any(Xorigin(zeroprob)>0) 
        negloglik = Inf;
        penalty = 2 * nparams;
        AIC = Inf;
    else
        if ~isempty(zeroprob)
            Phat = Phat(Phat>=1e-7);
        end
        negloglik = - 2 * sum(sum(Phat .* log(Phat)));
        penalty = 2 * nparams;
        AIC = negloglik + penalty;
    end
   

end

