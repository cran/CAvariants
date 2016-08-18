\name{sonscabasic}
\alias{sonscabasic}
\title{
Singly or one-way ordered non symmetrical correspondence analysis
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. 
This function is used in the main function \code{CAvariants} when the input parameter is \code{catype = "SONSCA"}.
It performs the hybrid decomposition of the numerator of the Goodman-Kruskal tau index and  implies an ordered (column) variable.
It compute the principal axes for rows and polynomial axes for columns, coordinates, weights of rows and columns, 
inertia (equal to the numerator of the tau index), the rank of the matrix.
It allows the analyst to decompose the inertia in column polynomial components. 
}
\usage{
sonscabasic(Xtable, mj)
}
\arguments{
\item{Xtable}{
The two-way contingency table.
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
mj<-c(1,2,3,4,5)
sonscabasic(asbestos,mj)
}
\keyword{ ~multivariate }
