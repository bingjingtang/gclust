\name{gclust.app}
\alias{gclust.app}
\concept{Graph Clustering}
\title{Graph Clustering Using SVT and NMF}
\description{Performs Clustering of Graphs using Singular Value Thresholding and Non-negative 
Factorization.}
\usage{
gclust.app(glist,r) 
}
\arguments{
\item{glist}{List of \code{igraph} Objects}
\item{r}{Maximum Number of Clusters Allowed}
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
}
