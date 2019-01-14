\name{nscabasic}
\alias{nscabasic}
\title{
Two-way non symmetrical correspondence analysis
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function is used in the main function \code{CAvariants} when the input parameter is \code{catype = "NSCA"}. 
It calculates the singular value decomposition of  the numerator of the Goodman-Kruskal tau index (index of predictability), 
computes principal axes, coordinates, weights of the rows and columns, total inertia (numerator of the tau index) and the rank of the matrix. 
}
\usage{
nscabasic(Xtable)
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
\examples{
data(asbestos)
nscabasic(asbestos)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate} % __ONLY ONE__ keyword per line
