\name{emerson.poly}
\alias{emerson.poly}
\title{
Orthogonal polynomials}
\description{
This function is called  from the functions
\code{docabasic, socabasic, sonscabasic} and \code{donscabasic}. 
It computes the orthogonal polynomials for the ordered categorical variables. 
The number of the polynomials is equal to the number of categories for that variable less one.
The function computes the polynomial transformation of the ordered categorical variable.}
\usage{
emerson.poly(mj, pj)
}
\arguments{
\item{mj}{
The ordered scores of the ordered variable. By default, \code{mj = NULL}, the natural scores (1, 2, ...) are computed.
}
\item{pj}{
The marginal relative frequencies of the ordered variable. 
}
}
\value{
Describe the value returned
\item{B}{the matrix of the orthogonal polynomials with the trivial polynomial removed.}
}
\references{
Beh EJ and Lombardo R 2014 Correspondence Analysis: Theory, Practice and New Strategies. Wiley.\cr
Emerson PL 1968  Numerical construction of orthogonal polynomials from a general recurrence formula. Biometrics, 24 (3), 695-701.\cr
Lombardo R Beh EJ 2016 Variants of Simple Correspondence Analysis. The R Journal, 8 (2), 167--184.\cr
Lombardo R Beh EJ and Kroonenberg PM 2016 Modelling Trends in Ordered Correspondence Analysis Using Orthogonal
Polynomials. Psychometrika, 81(2), 325-349.
}
\author{
Rosaria Lombardo and Eric J Beh
}
\note{
Note that the sum of the marginal relative frequencies of the ordered variables must be one.
}
%\examples{
%emerson.poly(c(1,2,3,4,5), as.vector(c(.1,.2,.3,.2,.2)))
%}
\keyword{nonparametric }% __ONLY ONE__ keyword per line
