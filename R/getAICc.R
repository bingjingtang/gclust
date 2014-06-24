getAICc <- function(gfit){
    Xorigin = gfit$Xorigin;
    Xmean = gfit$W %*% gfit$H %*% diag(colSums(Xorigin))
    Xmean[Xmean<1e-12] = 0
    Phat = gfit$W %*% gfit$H
    Phat[Phat < 1e-12] = 0
    if(ncol(gfit$W) > 1) { 
        nparams = colSums((1.0*(gfit$W > 0) ));nparams
        nparams = nparams %*% diag(1/colSums(Xorigin %*% t(gfit$H)));nparams
        nparams = sum(nparams);nparams
    } else {
        nparams = colSums((1.0*(gfit$W > 0) ));nparams
        nparams = nparams/colSums(Xorigin %*% t(gfit$H));nparams
        nparams = sum(nparams);nparams        
    }
    #     nparams = sum((1.0*(Phat > 0) ) %*% diag(1/colSums(Xorigin)));nparams
    zeroprob = which(Phat<1e-12)
    if(any(Xorigin[zeroprob]> 1e-12)) {
        retval = Inf
    } else {
        if(length(zeroprob) > 0) {
#            Xscaled = Xscaled[-zeroprob]
            Phat = Phat[-zeroprob]
        } 
        retval = 2*(1)*(nparams) - 2*sum(Phat*log(Phat))
    }
    data.frame(nclust=ncol(gfit$W),  negloglikpart=-2*sum(Phat*log(Phat)),  parampart=2*nparams*(1),   AIC=retval) 
}
