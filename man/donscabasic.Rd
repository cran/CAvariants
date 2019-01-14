\name{donscabasic}
\alias{donscabasic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Doubly, or two-way ordered, non symmetrical correspondence analysis: for two ordered variables}
\description{
This function is used in the main function \code{CAvariants} when the input parameter is \code{catype = "DONSCA"}. 
It performs the bivariate moment decomposition of the numerator of the Goodman-Kruskal tau index  for a contingency table consisting of two ordered variables.
It computes the polynomial axes, coordinates, weights of the rows and columns, total inertia (equal to the numerator of the  tau index) and the rank of the matrix. 
It also decomposes the inertia into row and column polynomial components. 
}
\usage{
donscabasic(Xtable, mi, mj)
}
\arguments{
\item{Xtable}{
The two-way contingency table.
}
\item{mi}{
The set of ordered row scores. By default, \code{mi = c(1:nrow(Xtable))} (natural scores).
}
\item{mj}{
The set of ordered  column scores. By default, \code{mj = c(1:ncol(Xtable))} (natural scores).
}
}
\references{
Beh EJ and Lombardo R 2014 Correspondence Analysis: Theory, Practice and New Strategies. Wiley.\cr
Lombardo R Beh EJ 2016 Variants of Simple Correspondence Analysis. The R Journal, 8 (2), 167--184.\cr
Lombardo R Beh EJ and Kroonenberg PM 2016 Modelling Trends in Ordered Correspondence Analysis Using Orthogonal
Polynomials. Psychometrika, 81(2), 325--349.

}
\author{
Rosaria Lombardo and Eric J. Beh
}
\examples{
data(asbestos)
mi <- c(1, 2, 3, 4) # natural scores for the rows
mj <- c(4.5, 14.5, 24.5, 34.5, 44.5) #midpoints for the columns
donscabasic(asbestos, mi, mj)
}
\keyword{ ~multivariate }
