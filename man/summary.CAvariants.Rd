\name{summary.CAvariants}
\alias{summary.CAvariants}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Summary of numerical results from CA variants 
}
\description{
This function prints a numerical summary of the results from any of the six variants of correspondence analysis. 
The input parameter is the name of the output of the main function  \code{CAvariants}. 
}
\usage{
\S3method{summary}{CAvariants}(object, printdims, digits, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{object}{
The output of the main function \code{CAvariants}.}
\item{printdims}{The number of dimensions that are used for summarising the numerical output of the analysis. By default, \code{printdims = 2}. 
the maximum number is equal to the rank of the table.}
\item{digits}{The number of decimal places used for displaying the numerical summaries of the analysis. 
By default, \code{digits = 3}.}
\item{...}{Further arguments passed to, or from, other functions.}
}
\value{
The value of output returned depends on the type of correspondence analysis that is performed.
\item{Inertias}{The inertia values, their percentage contribution to the total inertia and 
the cumulative percent inertias for the row and column variables.}
\item{Generalised correlation matrix}{The matrix of generalised correlations when performing 
an ordered correspondence analysis, \code{catype = "DOCA"}, \code{catype = "DONSCA"}, \code{catype = "SOCA"} or \code{catype = "SONSCA"}.}
\item{Row principal coordinates}{The row principal  coordinates when \code{catype = "CA"}, or \code{catype = "NSCA"}.}
\item{Column principal coordinates}{The column principal coordinates when \code{catype = "CA"}, or \code{catype = "NSCA"}.}
\item{Row standard coordinates }{The row standard coordinates when \code{catype = "CA"}, or \code{catype = "NSCA"}.}
\item{Column standard coordinates}{The column standard coordinates when \code{catype = "CA"}, or \code{catype = "NSCA"}.}
\item{Row principal polynomial coordinates}{The row principal polynomial coordinates 
when \code{catype = "DOCA"}, \code{catype = "DONSCA"}, \code{catype = "SOCA"}, or \code{catype = "SONSCA"}.}
\item{Column principal polynomial coordinates}{The column principal coordinates when \code{catype = "DOCA"}, or \code{catype = "DONSCA"}.}
\item{Row standard polynomial coordinates }{The row standard polynomial coordinates when \code{catype} is \code{"DOCA"} or \code{"DONSCA"}.}
\item{Column standard polynomial coordinates}{The column standard polynomial coordinates when \code{catype = "DOCA"}, \code{catype = "DONSCA"}, 
 \code{catype = "SOCA"}, or \code{catype = "SONSCA"}.}
\item{Total inertia}{The total inertia. For example, for non symmetrical correspondence analysis the numerator of the Goodman-Kruskal tau index, 
its C-statistic and p-value are returned.}
\item{Polynomial components}{The polynomial components of the total inertia and their p-values. 
The total inertia of the column space is partitioned to identify polynomial components. 
when \code{catype = "SOCA"} or \code{catype = "SONSCA"}. When \code{catype = "DOCA"} or 
\code{catype = "DONSCA"}, the total inertia of both the row and column space is partitioned 
to give the polynomial components.}
\item{Inner product}{The inner product of the biplot coordinates for the two-dimensional plot.}
}
\references{
Beh EJ and Lombardo R 2014 Correspondence Analysis: Theory, Practice and New Strategies. Wiley.\cr
Lombardo R Beh EJ 2016 Variants of Simple Correspondence Analysis. The R Journal, 8 (2), 167--184.\cr
Lombardo R Beh EJ and Kroonenberg PM 2016 Modelling Trends in Ordered Correspondence Analysis Using Orthogonal
Polynomials. Psychometrika, 81(2), 325--349.
}
\author{
Rosaria Lombardo and Eric J. Beh}
\examples{
asbestos <- matrix(c(310, 36, 0, 0, 212, 158, 9, 0, 21, 35, 17, 4, 25, 102,  
49, 18, 7, 35, 51, 28), 4, 5, dimnames = list(c("none","grade1", "grade2", "grade3"), 
c("0-9", "10-19", "20-29", "30-39", "40")))
risasbestos <- CAvariants(asbestos, catype = "DOCA", firstaxis = 1, lastaxis = 2) 
summary(risasbestos)
}
\keyword{summary}% __ONLY ONE__ keyword per line
