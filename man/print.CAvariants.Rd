\name{print.CAvariants}
\alias{print.CAvariants}
\title{
Main printing function for numerical summaries
}
\description{
This function prints the numerical output for any of the six variants of correspondence analysis called by \code{catype}. \cr
The input parameter is the name of the output of the main function  \code{CAvariants}. 
}
\usage{
\S3method{print}{CAvariants}(x, printdims = 2, ellcomp = TRUE, digits = 3,...) 
%print.CAvariants(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{The name of the output object from the main function \code{CAvariants}.}
\item{printdims}{The number of dimensions that are used for summarising the numerical output of the analysis. By default, \code{printdims = 2}. 
the maximum number is equal to the rank of the table.}
\item{ellcomp}{This parameter specifies whether the characteristics of the confidence ellipses (eccentricity, semi-axis, area, p-values) 
are to be computed. By default, \code{ellcomp =  TRUE}.  
}
\item{digits}{The number of decimal places used for displaying the numerical summaries of the analysis. 
By default, \code{digits = 3}.}
\item{...}{Further arguments passed to, or from, other functions.}
}
\details{
This function uses another function (called \code{printwithaxes}) for specifying the number of
columns of a matrix to print. 
}
\value{The output returned depends on the type of correspondence analysis that is performed
\item{Xtable}{The two-way contingency table.}
\item{Row weights: Imass}{The row weight matrix. These weights depend on the type of analysis that is performed.}
\item{Column weights: Jmass}{The column weight matrix. These weights are equal to the column marginal relative frequencies
for all types of analysis performed.}
\item{Total inertia}{The total inertia of the analysis performed. For example, 
for variants of non symmetrical correspondence analysis, the output produced includes the 
numerator of the Goodman-Kruskal tau index, its C-statistic and  p-value.}
\item{Inertias}{The inertia values, their percentage contribution to the total inertia and 
the cumulative percent inertias for the row and column variables.}
\item{Generalised correlation matrix}{The matrix of generalised correlations when performing 
an ordered correspondence analysis, \code{catype} must be \code{"DOCA"}, \code{"DONSCA"}, \code{"SOCA"} or \code{"SONSCA"}.}
\item{Row principal coordinates}{The row principal coordinates when \code{catype = "CA"} or \code{catype = "NSCA"}.}
\item{Column principal coordinates}{The column principal coordinates when \code{catype = "CA"} or \code{catype = "NSCA"}.}
\item{Row standard coordinates }{The row standard coordinates when \code{catype = "CA"} or \code{catype = "NSCA"}.}
\item{Column standard coordinates}{The column standard coordinates when \code{catype = "CA"} or \code{catype = "NSCA"}.}
\item{Row principal polynomial coordinates}{The row principal polynomial coordinates  when performing an ordered correspondence analysis.}
\item{Column principal polynomial coordinates}{The column principal coordinates  when performing a doubly  ordered correspondence analysis.}
\item{Row standard polynomial coordinates }{The row standard polynomial coordinates, when performing an ordered variant of correspondence analysis.}
\item{Column standard polynomial coordinates}{The column standard polynomial coordinates, when performing an ordered variant of correspondence analysis.}
\item{Row distances from the origin of the plot}{The squared Euclidean distance of the row categories from the origin of the plot.}
\item{Column distances from the origin of the plot}{The squared Euclidean distance of the column categories from the origin of the plot.}
\item{Polynomial components}{The polynomial components of the total inertia and their p-values. 
The total inertia of the column space is partitioned to identify polynomial components. 
when \code{catype = "SOCA"} or \code{catype = "SONSCA"}. When \code{catype = "DOCA"} or 
\code{catype = "DONSCA"},  the total inertia of both the row and column space is partitioned 
to give the polynomial components.}
%\item{Polynomial axes}{The column polynomial axes.}
%\item{Principal axes}{The row principal axes.}
\item{Inner product}{The inner product of the biplot coordinates for the two-dimensional plot.}
%\item{ellprint}{The logical parameter used for displaying the numerical characteristics of confidence ellipses. 
%By default, \code{ell = TRUE}.}
\item{eccentricity}{Value of ellipse eccentricity, the distance between its center and either of its two foci, It can be thought of as a measure of how much 
the conic section deviates from being circular.}
\item{HL Axis 1}{Value of ellipse semi-axis 1 for each row and column points.}
\item{HL Axis 2}{Value of ellipse semi-axis 2 for each row and column points.}
\item{Area}{Ellipse area  for each row and column points.}
\item{pvalcol}{P-value for each row and column points.}
}
\references{
Beh EJ and Lombardo R 2014 Correspondence Analysis: Theory, Practice and New Strategies. Wiley.\cr
Lombardo R Beh EJ 2016 Variants of Simple Correspondence Analysis. The R Journal, 8 (2), 167--184.
}
\author{
Rosaria Lombardo and Eric J. Beh}
\examples{
data(asbestos)
resasbestos <- CAvariants(asbestos, catype = "DOCA", firstaxis = 1, lastaxis = 2) 
print(resasbestos)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
