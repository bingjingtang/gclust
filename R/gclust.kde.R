reembed <- function(g, dmax, pamkout){
    X.list <- list()
    for(i in 1:pamkout$nc){
        idx.i <- which(pamkout$pamobject$clustering == i)
        gi <- g[idx.i,idx.i]
        Xhat.i <- adjacency.spectral.embedding(graph.adjacency(gi), dmax)$X
        ## eval <- sqrt(colSums(Xhat.i^2))
        ## dhat <- dimSelect(eval)
        X.list[[i]] <- Xhat.i
    }
    return(X.list)
}

rect.dist <- function(X,Y){
    n <- nrow(X)
    m <- nrow(Y)
    tmp1 <- X%*%t(Y)
    tmp2 <- outer(rep(1, n), rowSums(Y^2))
    tmp3 <- outer(rowSums(X^2), rep(1,m))
    
    D <- tmp2 - 2*tmp1 + tmp3
    return(D)
}

kernel.stat <- function(X,Y,sigma=0.2){
    n <- nrow(X)
    m <- nrow(Y)
    
    tmpXX <- sum(exp(-(as.matrix(dist(X))^2)/(2*sigma^2))) - n
    tmpYY <- sum(exp(-(as.matrix(dist(Y))^2)/(2*sigma^2))) - m
    tmpXY <- sum(exp(-(rect.dist(X,Y))/(2*sigma^2)))
    
    tmp <- tmpXX/(n*(n-1)) + tmpYY/(m*(m-1)) - 2*tmpXY/(m*n)
    
    return((m+n)*tmp)
}

find.transform <- function(X,Y){
    u <- apply(X,2,median)
    v <- apply(Y,2,median)
    if(ncol(X) == 1){
        T <- sign(u/v)
    }
    else{
        T <- diag(sign(u/v))
    }
    return(T)
}

gclust.kde <- function(X.list, sigma){
    K <- length(X.list)
    S <- matrix(0, K, K)
    for(i in 1:(K-1)){
        for(j in (i+1):K){
            Tij <- find.transform(X.list[[i]], X.list[[j]])
            S[i,j] <- kernel.stat(X.list[[i]] %*% Tij, X.list[[j]], sigma)
            S[j,i] <- S[i,j]
        }
    }
    return(S)
}
