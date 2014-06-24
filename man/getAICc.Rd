\name{getAICc}
\alias{getAICc}
\concept{Compute AIC based on a Poisson Approximation from \code{gclust}}
\title{Compute AIC based on a Poisson Approximation from \code{gclust} adjusting for finite sample}
\description{
Compute and Extract AIC information from \code{gclust} using a Poisson approximation,
where the penality term is adjusted for small sample cases.  Typically, 
\code{sbm.game} with small edge probability yields a such situation.  
}
\usage{ getAIC(x,rho) }
\arguments{
\item{x}{The output from a call to \code{gclust}}
}
\value{
Return a list vector containing the 4 numeric values for 
(1) the number of clulsters, 
(2) the negative log likelihood part of AIC, 
(3) the penality part of AIC, and 
(4) the value of AIC, in that order.
}

\author{ Nam Lee \email{nhlee@jhu.edu}}

\examples{
}
