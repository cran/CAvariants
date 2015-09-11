\name{summary.CAvariants}
\alias{summary.CAvariants}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Summary of results of CA variants analysis
}
\description{
This function prints summary results of any of the specified six variants of correspondence analysis. It is called internally from 
\code{CAvariants}. 
The input parameter is the name of the output of the main function  \code{CAvariants}. 
}
\usage{
\S3method{summary}{CAvariants}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{object}{
The output of the main function \code{CAvariants}.}
\item{...}{Further arguments passed to or from other methods.}
}
\value{
The value of output returned depends on the kind of correspondence analysis performed.
\item{Inertias}{The inertia value and its percentage, in the row and column space for each principal or polynomial axis when \code{catype} is \code{DOCA}, \code{DONSCA},
\code{SOCA} or \code{SONSCA}. When \code{catype} is \code{CA} or \code{NSCA},
the associated inertia values in the row and column spaces are the same.}
\item{Generalized correlation matrix}{The generalized correlation matrix when \code{catype} is \code{DOCA}, \code{DONSCA}, 
 \code{SOCA} or \code{SONSCA}.}
\item{Row principal coordinates}{The row principal  coordinates when \code{catype} is \code{CA} or \code{NSCA}.}
\item{Column principal coordinates}{The column principal coordinates when \code{catype} is \code{CA} or \code{NSCA}.}
\item{Row standard coordinates }{The row standard coordinates when \code{catype} is \code{CA} or \code{NSCA}.}
\item{Column standard coordinates}{The column standard coordinates when \code{catype} is \code{CA} or \code{NSCA}.}
\item{Row principal polynomial coordinates}{The row principal polynomial coordinates 
when \code{catype} is \code{DOCA}, \code{DONSCA}, \code{SOCA} or \code{SONSCA}.}
\item{Column principal plynomial coordinates}{The column principal coordinates when \code{catype} is \code{DOCA},\code{DONSCA}.}
\item{Row standard polynomial coordinates }{The row standard polynomial coordinates when \code{catype} is \code{DOCA} or \code{DONSCA}.}
\item{Column standard polynomial coordinates}{The column standard polynomial coordinates when \code{catype} is \code{DOCA}, \code{DONSCA}, 
 \code{SOCA}, \code{SONSCA}.}
\item{Total inertia}{The total inertia. For example in case of non symmetrical correspondence analysis, it gives
the numerator of the Goodman-Kruskal tau index, the associated C-statistic and its p-value.}
\item{Polynomial components}{The polynomial components of inertia and their p-values. The inertia of the column space is partitioned in terms
of polynomial components when \code{catype} is \code{SOCA} or \code{SONSCA}.
The inertia of the row and column space is partitioned in terms
of polynomial components when \code{catype} is \code{DOCA} or \code{DONSCA}}
\item{Inner product}{The inner product of coordinates.}
}
\references{
Beh EJ and Lombardo R 2014 Correspondence analysis, Theory, Practice and New Strategies. John Wiley & Sons.
}
\author{
Rosaria Lombardo and Eric J. Beh}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~
%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
asbestos<-matrix(c(310, 36, 0, 0, 212, 158, 9, 0, 21, 35, 17, 4, 25,102,  
49, 18, 7, 35, 51, 28), 4, 5, dimnames = list(c("none","grade1", "grade2", "grade3"), 
c("0-9", "10-19", "20-29", "30-39", "40")))
risasbestos<-CAvariants(asbestos,catype = "DOCA",firstaxis=1,lastaxis=2) 
summary.CAvariants(risasbestos)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
\keyword{ summary }% __ONLY ONE__ keyword per line
