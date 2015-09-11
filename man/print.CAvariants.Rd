\name{print.CAvariants}
\alias{print.CAvariants}
\title{
Main Printing function 
}
\description{
This function prints results of any of the specified six variants of correspondence analysis. It is called internally from 
\code{CAvariants}. 
The input parameter is the name of the output of the main function  \code{CAvariants}. 
}
\usage{
\S3method{print}{CAvariants}(x,printdims,...) 
%print.CAvariants(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{
The output of the main function \code{CAvariants}.}
\item{printdims}{ The dimension number in printing.
}
\item{...}{Further arguments passed to or from other methods.}
}
\details{
This function uses another function (called \code{printwithaxes}) for specifying the number of matrix dimensions to print. 
}
\value{
The value of output returned depends on the kind of correspondence analysis performed.
\item{DataTable }{The two-way contingency table.}
\item{Row weights }{The row weight matrix.}
\item{Column weights }{The column weight matrix.}
\item{Inertias}{The inertia value and its percentage, in the row and column space for each principal or polynomial axis when \code{catype} is \code{DOCA}, \code{DONSCA},
\code{SOCA} or \code{SONSCA}. When \code{catype} is \code{CA} or \code{NSCA},
the associated inertia values in the row and column spaces are the same.}
\item{Generalized correlation matrix}{The generalized correlation matrix when \code{catype} is \code{DOCA}, \code{DONSCA}, 
 \code{SOCA} or \code{SONSCA}.}
\item{Diagonal elements of the squared generalized correlation matrix Z'Z}{The diagonal elements of the squared generalized correlation matrix when \code{catype} is \code{DOCA}, \code{DONSCA}, 
 \code{SOCA} or \code{SONSCA}.}
\item{Diagonal elements of Z'Z for the constant (n-1)*(I-1)}{The column polynomial sources of inertia when \code{catype} is \code{DOCA}, \code{DONSCA}, 
 \code{SOCA} or \code{SONSCA}.}
\item{Diagonal elements of the squared generalized correlation matrix ZZ'}{The diagonal elements of the squared generalized correlation matrix  when \code{catype} is \code{DOCA}, \code{DONSCA}, 
 \code{SOCA} or \code{SONSCA}.}
\item{Diagonal elements of ZZ' for the constant (n-1)*(I-1)}{The row polynomial sources of inertia when \code{catype} is \code{DOCA}, \code{DONSCA}, 
 \code{SOCA} or \code{SONSCA}.}
\item{Row principal coordinates}{The row principal  coordinates when \code{catype} is \code{CA} or \code{NSCA}.}
\item{Column principal coordinates}{The column principal coordinates when \code{catype} is \code{CA} or \code{NSCA}.}
\item{Row standard coordinates }{The row standard coordinates when \code{catype} is \code{CA} or \code{NSCA}.}
\item{Column standard coordinates}{The column standard coordinates when \code{catype} is \code{CA} or \code{NSCA}.}
\item{Row principal polynomial coordinates}{The row principal polynomial coordinates when \code{catype} is \code{DOCA}, \code{DONSCA}, \code{SOCA} or \code{SONSCA}.}
\item{Column principal polynomial coordinates}{The column principal coordinates when \code{catype} is \code{DOCA},\code{DONSCA}.}
\item{Column distances from the origin of the plot}{The distance of the column categories from the origin of the plot.}
\item{Row distances from the origin of the plot}{The distance of the row categories from the origin of the plot.}
\item{Row standard polynomial coordinates }{The row standard polynomial coordinates when \code{catype} is \code{DOCA} or \code{DONSCA}.}
\item{Column standard polynomial coordinates}{The column standard polynomial coordinates when \code{catype} is \code{DOCA}, \code{DONSCA}, 
 \code{SOCA}, \code{SONSCA}.}
\item{Total inertia}{The total inertia. For example in case of non symmetrical correspondence analysis and its variants, it gives
the numerator of the Goodman-Kruskal tau index, the associated C-statistic and its p-value.}
\item{Polynomial axes}{The column polynomial axes.}
\item{Principal axes}{The row principal axes.}
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
print.CAvariants(risasbestos)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~multivariate }
%\keyword{multivariate}% __ONLY ONE__ keyword per line
