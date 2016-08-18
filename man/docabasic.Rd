\name{docabasic}
\alias{docabasic}
\title{
Doubly or two-way ordered correspondence analysis, it implies two ordered variables
}
\description{
This function is used in the main function \code{CAvariants} when the input parameter is \code{catype = "DOCA"}. 
It performs the bivariate moment decomposition of the  Pearson's ratio, 
computes polynomial axes, coordinates, weights of rows and columns, total inertia (equal to the Pearson's index), the rank of the matrix. 
It allows one to decompose the inertia into row and column polynomial components. 
}
\usage{
docabasic(Xtable, mi, mj)
}
\arguments{
\item{Xtable}{
The two-way contingency table.
}
\item{mi}{
The set of ordered row scores.
}
\item{mj}{
The set of ordered column scores.
}
}
\references{
Beh EJ and Lombardo R 2014 Correspondence Analysis: Theory, Practice and New Strategies. John Wiley & Sons.
}
\author{
Rosaria Lombardo and Eric J. Beh
}
\note{
This function belongs to the \code{R} object class called  \code{cabasicresults}.
}
\examples{
data(asbestos)
mi <- c(1,2,3,4)
mj <- c(1,2,3,4,5)
docabasic(asbestos, mi, mj)
}
\keyword{ ~multivariate }
