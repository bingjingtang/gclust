\name{gclust}
\alias{gclust}
\concept{Graph Clustering}
\title{Graph Clustering Using SVT and NMF}
\description{Performs Clustering of Graphs using Singular Value Thresholding and Non-negative 
Factorization.}
\usage{
gclust(glist,r,maxsvt=10) 
}
\arguments{
\item{glist}{List of \code{igraph} Objects}
\item{r}{Maximum Number of Clusters Allowed}
\item{maxsvt}{Maximum Number of Singular Value Thresholding to be done for estimating the mean matrix}
}
\value{
\item{nmf}{An NMF object}
\item{W}{Basis Graphs} 
\item{H}{Probability Vector for Cluster Weights}
\item{membership}{Integral Vector Specifying Cluster Membership Labels}
\item{AIC}{AIC Value based on an NMF factorizable Poisson Model (for a timeseries of graphs)}
}

\author{ Nam Lee \email{nhlee@jhu.edu}}

\examples{
###### sample usage  
n <- 100
k <- 2
m <- n/k
mvec <- rep(m,k)
cvec <- seq(0,0.5,by=0.1)
K <- length(cvec) # 6

## generate 6 basis graphs
set.seed(123456)
B <- lapply(1:K, function(x) (.5-cvec[x])*matrix(1,k,k)+(2*cvec[x])*diag(k))
g <- lapply(1:K, function(x) sbm.game(n,B[[x]],mvec))

## clone these 6 basis graph 100 times
R <- 100
Y <- sample(K,R,rep=T)
glist <- sapply(1:R, function(x) g[Y[[x]]])


## quick test run
#out = gclust(glist,6,10)
#length(unique(out$memb))
#adjustedRandIndex(Y,out$memb)

maxRank = 10
maxsvt = 10 
\dontrun{require(doMC)
registerDoMC(maxRank)
}
