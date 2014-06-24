getFICc <- function(gfit)  {
    X = gfit$Xorigin
    W = gfit$W
    H = gfit$H
    innDim = dim(W)[2]
    svdofX = svd(X)
    Uhat = svdofX$u[,1:innDim]
    Dhat = svdofX$d[1:innDim]
    Vhat = svdofX$v[,1:innDim]
    Xhat = W %*% H 
    if(innDim == 1) {
        Uhat = cbind(Uhat)
        Vhat = cbind(Vhat)
        Dhat = matrix(Dhat,1,1)
        dW = W - X %*% t(Xhat) %*% Uhat %*% matrix(1/Dhat^2,1,1) %*% t(Uhat) %*% W
        dH = t(H) - t(X) %*% Xhat %*% Vhat %*% matrix(1/Dhat^2,1,1) %*% t(Vhat) %*% t(H)
    } else {
        dW = W - X %*% t(Xhat) %*% Uhat %*% diag(1/Dhat^2) %*% t(Uhat) %*% W
        dH = t(H) - t(X) %*% Xhat %*% Vhat %*% diag(1/Dhat^2) %*% t(Vhat) %*% t(H)
    }
    retval = (c(norm(X-Xhat,type='F'),norm(dW,type='F'),norm(dH,type='F')))
    data.frame(nclust=innDim,  negloglikpart=retval[1],  parampart=sum(retval[2:3]),   FIC=sum(retval))
}