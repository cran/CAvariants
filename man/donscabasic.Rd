\name{donscabasic}
\alias{donscabasic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Doubly or two-way ordered non symmetrical correspondence analysis}
\description{
This function is used in the main function \code{CAvariants} when the input parameter is \code{catype = "DONSCA"}. 
It performs the bivariate moment decomposition of the numerator of the Goodman-Kruskal tau index  and implies two ordered variables.
It computes polynomial axes, coordinates, weights of rows and columns, total inertia (equal to the numerator of the  tau index) and the rank of the matrix. 
It allows the analyst to decompose the inertia into row and column polynomial components. 
}
\usage{
donscabasic(Xtable, mi, mj)
}
\arguments{
\item{Xtable}{
The two-way contingency table.
}
\item{mi}{
The set of ordered row scores.
}
\item{mj}{
The set of ordered  column scores.
}
}
\references{
Beh EJ and Lombardo R 2014 Correspondence Analysis: Theory, Practice and New Strategies. John Wiley & Sons.\cr
Lombardo R Beh EJ 2016 Variants of Simple Correspondence Analysis. The R Journal, 8 (2), 167--184.\cr
Lombardo R Beh EJ and Kroonenberg PM 2016 Modelling Trends in Ordered Correspondence Analysis Using Orthogonal
Polynomials. Psychometrika, 81(2), 325--349.

}
\author{
Rosaria Lombardo and Eric J. Beh
}
\note{
This function belongs to the \code{R} object class called \code{cabasicresults}.
}
\examples{
data(asbestos)
mi <- c(1,2,3,4)
mj <- c(1,2,3,4,5)
donscabasic(asbestos, mi, mj)
}
\keyword{ ~multivariate }
