\name{getFICc}
\alias{getFICc}
\concept{Compute FIC based on an fixed point error formula, from \code{gclust}}
\title{Compute FIC based on an fixed point error formula, from \code{gclust} adjusting for finite sample}
\description{
Compute and Extract FIC information from \code{gclust} using a Poisson approximation,
where the penality term is adjusted for small sample cases.  Typically, 
\code{sbm.game} with small edge probability yields a such situation.  
}
\usage{ getFIC(x) }
\arguments{
\item{x}{The output from a call to \code{gclust}}
}
\value{
Return a list vector containing the 4 numeric values for 
(1) the number of clulsters, 
(2) the negative log likelihood part of FIC, 
(3) the penality part of FIC, and 
(4) the value of FIC, in that order.
}

\author{ Nam Lee \email{nhlee@jhu.edu}}

\examples{
}
