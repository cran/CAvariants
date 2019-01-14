\name{sonscabasic}
\alias{sonscabasic}
\title{
Singly, or one-way, ordered non symmetrical correspondence analysis: for an ordered column predictor variable
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. 
This function is used by the main function \code{CAvariants} when the input parameter is \code{catype = "SONSCA"}.
It performs the hybrid decomposition of the numerator of the Goodman-Kruskal tau index and  implies an ordered (column) variable.
It calculates the principal axes for the rows and polynomial axes for the columns, coordinates.
It also calculates the row and column weights, inertia (based on the numerator of the tau index) and the rank of the matrix.
It decomposes the inertia into column polynomial components. 
}
\usage{
sonscabasic(Xtable, mj)
}
\arguments{
\item{Xtable}{
The two-way contingency table.
}
\item{mj}{
The set of ordered column scores. By default, \code{mj = c(1:ncol(Xtable)) (natural scores).}}
}
\references{
Beh EJ and Lombardo R 2014 Correspondence Analysis: Theory, Practice and New Strategies. Wiley.\cr
Lombardo R Beh EJ 2016 Variants of Simple Correspondence Analysis. The R Journal, 8 (2), 167--184.\cr
Lombardo R Beh EJ and Kroonenberg PM 2016 Modelling Trends in Ordered Correspondence Analysis Using Orthogonal
Polynomials. Psychometrika, 81(2), 325-349.
}
\author{
Rosaria Lombardo and Eric J. Beh
}
\examples{
data(asbestos)
mj<-c(1, 2, 3, 4, 5)
sonscabasic(asbestos, mj)
}
\keyword{multivariate}
\keyword{orthogonal polynomials}
