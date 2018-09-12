\name{summary.CAvariants}
\alias{summary.CAvariants}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Summary of results of CA variants analysis
}
\description{
This function prints summary results of any of the specified six variants of correspondence analysis. 
The input parameter is the name of the output of the main function  \code{CAvariants}. 
}
\usage{
\S3method{summary}{CAvariants}(object,printdims,digits,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{object}{
The output of the main function \code{CAvariants}.}
\item{printdims}{The dimension number in printing. By default, \code{printdims = 3}.}
\item{digits}{The minimum number of significant digits to be printed in values. By default, \code{digits = 3}.}
\item{...}{Further arguments passed to or from other methods.}
}
\value{
The value of output returned depends on the kind of correspondence analysis performed.
\item{Inertias}{The inertia value and its percentage, in the row and column space for each principal or polynomial axis when \code{catype} is \code{"DOCA"}, \code{"DONSCA"},
\code{"SOCA"} or \code{"SONSCA"}. When \code{catype} is \code{"CA"} or \code{"NSCA"},
the associated inertia values in the row and column spaces are the same.}
\item{Generalized correlation matrix}{The generalized correlation matrix when \code{catype} is \code{"DOCA"}, \code{"DONSCA"}, 
 \code{"SOCA"} or \code{"SONSCA"}.}
\item{Row principal coordinates}{The row principal  coordinates when \code{catype} is \code{"CA"} or \code{"NSCA"}.}
\item{Column principal coordinates}{The column principal coordinates when \code{catype} is \code{"CA"} or \code{"NSCA"}.}
\item{Row standard coordinates }{The row standard coordinates when \code{catype} is \code{"CA"} or \code{"NSCA"}.}
\item{Column standard coordinates}{The column standard coordinates when \code{catype} is \code{"CA"} or \code{"NSCA"}.}
\item{Row principal polynomial coordinates}{The row principal polynomial coordinates 
when \code{catype} is \code{"DOCA"}, \code{"DONSCA"}, \code{"SOCA"} or \code{"SONSCA"}.}
\item{Column principal polynomial coordinates}{The column principal coordinates when \code{catype} is \code{"DOCA"},\code{"DONSCA"}.}
\item{Row standard polynomial coordinates }{The row standard polynomial coordinates when \code{catype} is \code{"DOCA"} or \code{"DONSCA"}.}
\item{Column standard polynomial coordinates}{The column standard polynomial coordinates when \code{catype} is \code{"DOCA"}, \code{"DONSCA"}, 
 \code{"SOCA"}, \code{"SONSCA"}.}
\item{Total inertia}{The total inertia. For example in case of non symmetrical correspondence analysis, it gives
the numerator of the Goodman-Kruskal tau index, the associated C-statistic and its p-value.}
\item{Polynomial components}{The polynomial components of inertia and their p-values. The inertia of the column space is partitioned in terms
of polynomial components when \code{catype} is \code{"SOCA"} or \code{"SONSCA"}.
The inertia of the row and column space is partitioned in terms
of polynomial components when \code{catype} is \code{"DOCA"} or \code{"DONSCA"}.}
\item{Inner product}{The inner product of coordinates.}
}
\references{
Beh EJ and Lombardo R 2014 Correspondence Analysis: Theory, Practice and New Strategies. John Wiley & Sons.\cr
Lombardo R Beh EJ 2016 Variants of Simple Correspondence Analysis. The R Journal, 8 (2), 167--184.\cr
Lombardo R Beh EJ and Kroonenberg PM 2016 Modelling Trends in Ordered Correspondence Analysis Using Orthogonal
Polynomials. Psychometrika, 81(2), 325--349.
}
\author{
Rosaria Lombardo and Eric J. Beh}
\examples{
asbestos<-matrix(c(310, 36, 0, 0, 212, 158, 9, 0, 21, 35, 17, 4, 25,102,  
49, 18, 7, 35, 51, 28), 4, 5, dimnames = list(c("none","grade1", "grade2", "grade3"), 
c("0-9", "10-19", "20-29", "30-39", "40")))
risasbestos<-CAvariants(asbestos, catype = "DOCA", firstaxis = 1, lastaxis = 2) 
summary(risasbestos)
}
\keyword{ summary }% __ONLY ONE__ keyword per line
