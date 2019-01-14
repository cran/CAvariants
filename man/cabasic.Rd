\name{cabasic}
\alias{cabasic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Classical two-way correspondence analysis
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function is used in the main function \code{CAvariants} when the input parameter is \code{catype = "CA"}. 
It performs the singular value decomposition of Pearson's ratio and
computes principal axes, coordinates, the weights of rows and columns, the total inertia (equal to Pearson's index) 
and the rank of the matrix. 
}
\usage{
cabasic(Xtable)
}
\arguments{
\item{Xtable}{
The two-way contingency table.
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
\note{
This function belongs to the \code{R} object class called  \code{cabasicresults}.
}
\examples{
data(asbestos)
cabasic(asbestos)
}
\keyword{ multivariate}
