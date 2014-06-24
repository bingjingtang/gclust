gclust.rsvt <- function(glist,r=1,maxsvt=10,nmfout=FALSE,maxit=10000, nmfmethod='lee',...) 
{
    # Model Based Graph Clutering Approach    
    if (is.list(glist)) {
        Time <- length(glist)
        if (class(glist[[1]])=="igraph") {
            glist <- lapply(glist,get.adjacency)
        }
        Xorigin <- sapply(glist,as.vector)
    } else {
        Xorigin = glist
    }
    
    Xraw = Xorigin 
    Xorigin = Xorigin %*% solve(diag(colSums(Xorigin)))
    mysvd = tryCatch(irlba(Xorigin,r,r,maxit=maxit),error=function(e) svd(Xorigin,r,r))
    U = mysvd$u[,1:r,drop=FALSE]
    V = mysvd$v[,1:r,drop=FALSE]
    if(r == 1) {
        S = matrix(mysvd$d[1:r],1,1)
    } else {
        S = diag(mysvd$d[1:r])
    }
    
    maxsvt = ifelse(maxsvt > 0,maxsvt,1)
    for(itr in 1:maxsvt) {
        X = U %*% S %*% t(V)
        X[X<0]=0
        mysvd = tryCatch(irlba(X,r,r,maxit=maxit),error=function(e) svd(X,r,r))
        UU = mysvd$u[,1:r,drop=FALSE]
        VV = mysvd$v[,1:r,drop=FALSE]
        if(r == 1) {
            SS = matrix(mysvd$d[1:r],1,1)
        } else {
            SS = diag(mysvd$d[1:r])
        }
        if(norm(UU-U) + norm(VV-V) + norm(SS-S) < 1e-12) {
            break
        }
        U = UU;
        V = VV;
        S = SS;
    }
    
    XX = X %*% solve(diag(colSums(X)))
    stash = which(rowSums(XX)==0)
    if(length(stash)>0){
        XX =  XX[-stash,]
    }
    mynmf = nmf(XX,rank=r,method=nmfmethod,...)
    WW = matrix(0,nrow=nrow(X),ncol=r)
    if(length(stash)>0){
        WW[-stash,] = basis(mynmf)
    } else {
        WW = basis(mynmf)
    }

    HH = coef(mynmf)
    if(r > 1) {
        HH = diag(colSums(WW)) %*% HH
        WW = round(WW %*% solve(diag(colSums(WW))),12)
        HH = round(HH %*% solve(diag(colSums(HH))),12)   
    } else {
        HH = colSums(WW) %*% HH
        WW = round(WW /colSums(WW),12)
        HH = round(HH %*% solve(diag(colSums(HH))),12)   
    }
    
    if(nmfout)
        return(list(nmf=mynmf,W=WW,H=HH, Xorigin=Xraw, Xmean=X %*% diag(colSums(Xraw))))
    else 
        return(list(nmf=NULL,W=WW,H=HH, Xorigin=Xraw, Xmean=X %*% diag(colSums(Xraw))))
}


gclust.app <- function(glist,r=1, nmfout=FALSE,maxit=10000,nmfmethod='lee',...) 
{
    # Model Based Graph Clutering Approach    
    if (is.list(glist)) {
        Time <- length(glist)
        if (class(glist[[1]])=="igraph") {
            glist <- lapply(glist,get.adjacency)
        }
        Xorigin <- sapply(glist,as.vector)
    } else {
        Xorigin = glist
    }
    
    Xraw = Xorigin 
    XX = Xorigin %*% solve(diag(colSums(Xorigin)))
    
    stash = which(rowSums(XX)==0)
    if(length(stash)>0){
        XX =  XX[-stash,]
    }
    mynmf = nmf(XX,rank=r,method=nmfmethod,...)
    WW = matrix(0,nrow=nrow(Xraw),ncol=r)
    if(length(stash)>0){
        WW[-stash,] = cbind(basis(mynmf))
    } else {
        WW = basis(mynmf)
    }
    HH = coef(mynmf)
    if(r > 1) {
        HH = diag(colSums(WW)) %*% HH
        WW = round(WW %*% solve(diag(colSums(WW))),12)
        HH = round(HH %*% solve(diag(colSums(HH))),12)   
    } else {
        HH = colSums(WW) %*% HH
        WW = round(WW /colSums(WW),12)
        HH = round(HH %*% solve(diag(colSums(HH))),12)   
    }
    
    if(nmfout)
        return(list(nmf=mynmf,W=WW,H=HH, Xorigin = Xraw, Xmean=Xraw %*% diag(colSums(Xraw))))
    else 
        return(list(nmf=NULL,W=WW,H=HH, Xorigin = Xraw, Xmean=Xraw %*% diag(colSums(Xraw))))
}

