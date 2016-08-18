\name{print.CAvariants}
\alias{print.CAvariants}
\title{
Main printing function 
}
\description{
This function prints the results of any of the specified six variants of correspondence analysis. \cr
It is called internally from  \code{CAvariants}. The input parameter is the name of the output of the main function  \code{CAvariants}. 
}
\usage{
\S3method{print}{CAvariants}(x, printdims = 2, ellprint = TRUE, Mell = min(nrow(x$Xtable), ncol(x$Xtable)) - 1, 
 alpha = 0.05, digits = 3,...) 
%print.CAvariants(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{The name of the output object, for example say \code{res}, used with the main function \code{CAvariants}.}
\item{printdims}{ The number of dimensions, \code{printdims}, that are used to generate the correspondence plot, or biplot, 
and for summarising the numerical output of the analysis. By default, \code{printdims  =  2}.}
\item{ellprint}{The flag parameter, \code{ellprint}, allows that the characteristics of the confidence ellipses (eccentricity, semi-axis, area, p-values)  are displayed. \cr
By default, \code{ellprint  =  TRUE}. }
\item{Mell}{The number of axes \code{Mell} used for the construction of the confidence ellipses. By default, 
it is equal to its maximum value, \code{Mell = min(nrow(x$Xtable), ncol(x$Xtable)) - 1}, i.e. the rank of the data matrix.}
\item{alpha}{The level of significance used for the construction of the elliptical regions, \code{alpha}. By default, \code{alpha = 0.05}.}
\item{digits}{The minimum number of decimal places, \code{digits}, used for displaying the numerical summaries of the analysis. 
By default, \code{digits = 3}.}
\item{...}{Further arguments passed to or from other methods.}
}
\details{
This function uses another function (called \code{printwithaxes}) for specifying the number of matrix dimensions to print. 
}
\value{The value of output returned depends on the kind of correspondence analysis performed
\item{Xtable}{The two-way contingency table.}
\item{Row weights: Imass}{The row weight matrix. These weights depend on the type of analysis performed.}
\item{Column weights: Jmass}{The column weight matrix. These weights are equal to the data column margins for all types of analysis performed.}
\item{Total inertia}{The total inertia of the analysis performed. For example, when considering the variants \cr
of non symmetrical correspondence analysis, the numerator of the Goodman-Kruskal tau index, \cr
the associated C-statistic and its p-value is produced.}
\item{Inertias}{The inertia values, their percentage contribution to the total inertia and the cumulative \cr
percent inertias of the row and column space \code{Inertias}. \cr
When performing an ordered correspondence analysis, this output summary describes both the row and column spaces \cr
for each principal or polynomial axis. When \code{catype} is \code{CA} or \code{NSCA}, \cr
the associated inertia values in the row and column spaces are identical.}
\item{Generalized correlation matrix}{The generalized correlation matrix when performing an ordered correspondence analysis, \cr
\code{catype} should be \code{"DOCA"}, \code{"DONSCA"},  \code{"SOCA"} or \code{"SONSCA"}.}
\item{Row principal coordinates}{The row principal  coordinates when \code{catype} is \code{"CA"} or \code{"NSCA"}.}
\item{Column principal coordinates}{The column principal coordinates when \code{catype} is \code{"CA"} or \code{"NSCA"}.}
\item{Row standard coordinates }{The row standard coordinates when \code{catype} is \code{"CA"} or \code{"NSCA"}.}
\item{Column standard coordinates}{The column standard coordinates when \code{catype} is \code{"CA"} or \code{"NSCA"}.}
\item{Row principal polynomial coordinates}{The row principal polynomial coordinates  when performing an ordered correspondence analysis.}
\item{Column principal polynomial coordinates}{The column principal coordinates  when performing a doubly  ordered correspondence analysis.}
\item{Row standard polynomial coordinates }{The row standard polynomial coordinates, i.e. row polynomial axes, when performing a doubly ordered correspondence analysis.}
\item{Column standard polynomial coordinates}{The column standard polynomial coordinates, i.e. column polynomial axes,  when performing an ordered correspondence analysis.}
\item{Row distances from the origin of the plot}{The Euclidean distance of the row categories from the origin of the plot.}
\item{Column distances from the origin of the plot}{The Euclidean distance of the column categories from the origin of the plot.}
\item{Polynomial components}{The polynomial components of the total inertia and their p-values. The total inertia of the column space \cr
is partitioned to identify polynomial components, when \code{catype} is \code{"SOCA"} or \code{"SONSCA"}. When \code{catype}  \cr
is \code{"DOCA"} or \code{"DONSCA"},  the total inertia of both the row and column space is partitioned to identify of polynomial components.}
%\item{Polynomial axes}{The column polynomial axes.}
%\item{Principal axes}{The row principal axes.}
\item{Inner product}{The inner product of the biplot coordinates (concerning the first two axes when \code{firstaxis = 1} and \code{lastaxis = 2}). }
\item{ellprint }{When the input flag parameter is  \code{ellprint = TRUE}, then the print includes the  eccentricity of the confidence ellipses,  
the semi-major axis length of the ellipse for each row and column point, \code{HL Axis 1}, the semi-minor axis length for the ellipse for each row and column point, 
\code{HL Axis 2}, the area of the ellipse for each row and column point, \code{Area} and the p-value for each row and column point, \code{P-value}, 
see also the parameter \code{ellcomp} of the function \code{CAvariants} for a detailed description of these parameters.
}
%\item{eccentricity}{Value of ellipse eccentricity, the distance between its center and either of its two foci, It can be thought of as a measure of how much 
%the conic section deviates from being circular.}
%\item{HL Axis 1}{Value of ellipse semi-axis 1 for each row and column points.}
%\item{HL Axis 2}{Value of ellipse semi-axis 2 for each row and column points.}
%\item{Area }{Ellipse area  for each row and column points.}
%\item{pvalcol }{P-value for each row and column points.}
}
\references{
Beh EJ and Lombardo R 2014 Correspondence Analysis: Theory, Practice and New Strategies. John Wiley & Sons.
}
\author{
Rosaria Lombardo and Eric J. Beh}
\examples{
data(asbestos)
risasbestos<-CAvariants(asbestos, catype = "DOCA", firstaxis = 1, lastaxis = 2) 
print.CAvariants(risasbestos)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
