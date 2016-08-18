\name{compstable.exe}
\alias{compstable.exe}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Polynomial component of inertia in row and column spaces
}
\description{
This function allows the analyst to compute the contribution of the polynomial components to the inertia (chi-squared or tau).
The ordered variable should be the both row and column variables that are transformed by polynomials.
The polynomial components are the row and column polynomial components. 
The given input matrix is the Z matrix of generalised correlations from the bivariate moment decomposition.
It is called by \code{CAvariants} when \code{catype="DOCA"} or \code{catype = "DONSCA"}.
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
\item{comps }{The matrix of the polynomial component of inertia.}
}
\references{
Beh EJ and Lombardo R 2014 Correspondence Analysis: Theory, Practice and New Strategies. John Wiley & Sons.
}
\author{
Rosaria Lombardo and Eric J. Beh
}
\note{
This function belongs to the \code{R} object class called  \code{cacorporateplus}.
}
\keyword{ multivariate }% __ONLY ONE__ keyword per line
