\name{compstable.exe}
\alias{compstable.exe}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Polynomial component of inertia for the row and column spaces
}
\description{
This function allows the analyst to compute the contribution of the polynomial components to the inertia 
(Pearson's chi-squared statistic or the Goodman-Kruskal tau index).
The ordered variable should be both the row and column variables that are transformed by the polynomials.
The polynomial components are the row and column polynomial components. 
The given input matrix is the Z matrix of generalised correlations from the bivariate moment decomposition.
It is called by \code{CAvariants} when \code{catype = "DOCA"} or \code{catype = "DONSCA"}.
}
\usage{
compstable.exe(Z)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{Z}{The matrix of generalised correlations between the polynomial axes.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
The value returned is the matrix
\item{comps }{The matrix of the polynomial components of the inertia.}
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
\keyword{ multivariate }% __ONLY ONE__ keyword per line
