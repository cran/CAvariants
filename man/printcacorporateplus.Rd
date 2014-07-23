\name{printcacorporateplus}
\alias{printcacorporateplus}
\title{
Main Printing function 
}
\description{
This function prints results of any of the specified six variants of correspondence analysis. It is called internally from 
\code{CAvariants}. 
The input parameter is the name of the output of the main function  \code{CAvariants}. 
}
\usage{
printcacorporateplus(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{
The output of the main function \code{CAvariants}.}
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
\item{Row principal polynomial coordinates}{The row principal polynomial coordinates 
when \code{catype} is \code{DOCA}, \code{DONSCA}, \code{SOCA} or \code{SONSCA}.}
\item{Column principal plynomial coordinates}{The column principal coordinates when \code{catype} is \code{DOCA},\code{DONSCA}.}
\item{Row standard polynomial coordinates }{The row standard polynomial coordinates when \code{catype} is \code{DOCA} or \code{DONSCA}.}
\item{Column standard polynomial coordinates}{The column standard polynomial coordinates when \code{catype} is \code{DOCA}, \code{DONSCA}, 
 \code{SOCA}, \code{SONSCA}.}
\item{Total inertia}{The total inertia. For example in case of non symmetrical correspondence analysis, it gives
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
%\examples{
%print.cacorporateplus(x)
%function (x) 
%{
%    d <- min(x@printdims, x@S@r)
%    axnames <- character(length = d)
%    for (i in 1:d) {
%        axnames[i] <- paste(" Axis", i)
%    }
%    cat("\n    RESULTS for Corporate Correspondence Analysis:\n")
%    cat("\n    Data Table:\n")
%    print(x@DataMatrix)
%    cat("\n    Row Weights: Imass:\n")
%    round(print(matrix(x@Imass, x@rows, x@rows, dimnames = list(x@rowlabels, 
%        x@rowlabels)), dig = 3))
%    cat("\n    Column Weights: Jmass:\n")
%    round(print(matrix(x@Jmass, x@cols, x@cols, dimnames = list(x@collabels, 
%        x@collabels)), dig = 3))
%    cat("Total inertia ", round(x@inertiasum, dig = 3), "\n\n")
%    cat("Inertias, percent inertias and cumulative percent inertias \n\n")
%    print(data.frame(x@inertias))
%    cat("\nRows in standard coordinates\n\n")
%    printwithaxes(data.frame(x@Rstdcoord[, 1:d], row.names = x@rowlabels), 
%        axnames)
%    cat("\nRows in principal coordinates\n\n")
%    printwithaxes(data.frame(x@Rprinccoord[, 1:d], row.names = x@rowlabels), 
%        axnames)
%    cat("\nColumns in principal coordinates\n\n")
%    printwithaxes(data.frame(x@Cprinccoord[, 1:d], row.names = x@collabels), 
%        axnames)
%    cat("\nColumns in standard coordinates\n\n")
%    printwithaxes(data.frame(x@Cstdcoord[, 1:d], row.names = x@collabels), 
%        axnames)
%    if ((x@catype == "NSCA") || (x@catype == "DONSCA") || (x@catype == 
%        "SONSCA")) {
%        cat("\n    Predictability Index for variants of Non symmetrical Correspondence Analysis:\n")
%        cat("\nTau Index predicting from column \n\n")
%        print(x@S@tau)
%        Cstatistic <- (sum(x@DataMatrix) - 1) * (nrow(x@DataMatrix) - 
%            1) * x@S@tau
%        pvalueC <- 1 - pchisq(Cstatistic, (nrow(x@DataMatrix) - 
%            1) * (ncol(x@DataMatrix) - 1))
%        cat("\n C-statistic", Cstatistic, "and p-value", pvalueC, 
%            "\n")
%    }
%    if ((x@catype == "DOCA") | (x@catype == "DONSCA")) {
%        cat("\n Column Polynomial axes \n")
%        printwithaxes(data.frame(x@S@Raxes[, 1:d], row.names = x@collabels), 
%            axnames)
%        cat("\n Row Polynomial axes \n")
%        printwithaxes(data.frame(x@S@Caxes[, 1:d], row.names = x@rowlabels), 
%            axnames)
%    }
%    if ((x@catype == "SOCA") | (x@catype == "SONSCA")) {
%        cat("\n Column Polynomial axes \n")
%        printwithaxes(data.frame(x@S@Raxes[, 1:d], row.names = x@collabels), 
%            axnames)
%        cat("\n Row Principal axes \n")
%        printwithaxes(data.frame(x@S@Caxes[, 1:d], row.names = x@rowlabels), 
%            axnames)
%    }
%    else {
%        cat("\n Column Principal axes \n")
%        printwithaxes(data.frame(x@S@Raxes[, 1:d], row.names = x@collabels), 
%            axnames)
%        cat("\n Row Principal axes \n")
%        printwithaxes(data.frame(x@S@Caxes[, 1:d], row.names = x@rowlabels), 
%            axnames)
%    }
%    if ((x@catype == "DOCA") || (x@catype == "SOCA") || (x@catype == 
%        "DONSCA") || (x@catype == "SONSCA")) {
%        cat("\n Polynomial Components of Inertia \n")
%        print(x@comps)
%    }
%    else {
%        cat("\n\n")
%    }
%    cat("\n Inner Product of coordinates (first two axes)   \n")
%    print(round(x@Trend, dig = 3))
%  }
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~multivariate }
%\keyword{multivariate}% __ONLY ONE__ keyword per line
