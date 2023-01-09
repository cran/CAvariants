\name{CAvariants}
\alias{CAvariants}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Six variants of correspondence analysis 
}
\description{
It performs \cr
1) simple correspondence analysis \cr
2) doubly ordered correspondence analysis \cr
3) singly ordered correspondence analysis  \cr
4) non symmetrical correspondence analysis \cr
5) doubly ordered non symmetrical correspondence analysis \cr
6) singly ordered non symmetrical correspondence analysis 
}
\usage{
CAvariants(Xtable, mj = NULL, mi = NULL, firstaxis = 1, lastaxis = 2,
catype = "CA",  M = min(nrow(Xtable), ncol(Xtable)) - 1, alpha = 0.05) 
}
\arguments{
\item{Xtable}{
The two-way contingency table.
}
\item{mi}{
The assigned ordered scores for the row categories. By default, \code{mi = NULL}, which gives consecutive integer valued (natural) scores.
}
\item{mj}{ 
The assigned ordered scores for the column categories,  By default, \code{mj = NULL}, which gives consecutive integer valued (natural) scores.
} 
%\item{printdims}{
%The parameter for specifying the dimension of the coordinate and distance matrices in print. By default  \code{printdims = 2}.
%}
\item{firstaxis}{
The horizontal polynomial, or principal, axis.  It is used for the construction of the Inner product table. By default \code{firstaxis = 1}.
}
\item{lastaxis}{
The vertical polynomial, or principal, axis. It is used  for the construction of the Inner product table. By default \code{lastaxis = 2}.
}
\item{catype}{
The input parameter for specifying what variant of correspondence analysis is to be performed. By default,  \code{catype = "CA"}. 
Other possible values are:  \code{catype = "SOCA"}, \code{catype = "DOCA"}, \code{catype = "NSCA"}, \code{catype = "SONSCA"}, \code{catype = "DONSCA"}.
}
%\item{ellcomp}{This input parameter ensures that the characteristics of the algebraic confidence ellipses are computed and stored. 
%When \code{ellcomp = TRUE} (which is by default), the output includes the characteristics of the ellipses. 
%The eccentricity of the confidence ellipses is summarised by the quantity \code{eccentricity}, which is the distance between the center and either 
%of its two foci,
% and can be thought of as a measure of how much the conic section deviates from being circular
%(when the regions are perfectly circular than the eccentricity is zero).
% The semi-major axis length of the ellipse for each row and column 
%point is given by \code{HL Axis 1} while \code{HL Axis 2} gives the semi-minor axis length of the points. 
%The area of the ellipse for each row and column category is given by \code{Area} while the p-value of each category  is given by \code{P-value}.}
\item{M}{
The number of axes used for determining the structure of the elliptical confidence regions. 
By default,  \code{M = min(nrow(Xtable), ncol(Xtable)) - 1}, i.e. the rank of the data matrix.
}
\item{alpha}{
The level of significance for the elliptical regions. By default, \code{alpha = 0.05}.
}
}
%\details{
%This function belongs to the object class called \code{cacorporateplus}  
%}
\value{
Description of the output returned
\item{Xtable}{The two-way contingency table.}
\item{rows}{The number of rows of the two-way contingency table.}
\item{cols}{The number of columns of the two-way contingency table.}
\item{r}{The rank of the two-way contingency table.}
\item{n}{The total number of observations of the two-way contingency table.}
\item{rowlabels}{The labels of the row variable.}
\item{collabels}{The labels of the column variable.}
\item{Rprinccoord}{The row principal  coordinates.  When the input parameter \code{catype} is \cr
\code{"DOCA"}, \code{"SOCA"}, \code{"DONSCA"} or \code{"SONSCA"}, they are row principal polynomial coordinates.}
\item{Cprinccoord}{The column principal  coordinates. When the input parameter \code{catype} is  \cr
\code{"DOCA"}, \code{"SOCA"}, \code{"DONSCA"} or \code{"SONSCA"}, they are column principal polynomial coordinates.}
\item{Rstdcoord }{The row standard  coordinates.  When the input parameter \code{catype} is \cr
\code{"DOCA"}, \code{"SOCA"}, \code{"DONSCA"} or \code{"SONSCA"}, they are row standard polynomial coordinates.}
\item{Cstdcoord }{The column standard  coordinates. When the input parameter  \code{catype} is  \cr
\code{"DOCA"}, \code{"SOCA"}, \code{"DONSCA"} or \code{"SONSCA"}, they are column standard polynomial coordinates.}
\item{tauden}{The denominator of the Goodman-Kruskal tau index is given when the input parameter  \code{catype}  is \code{"NSCA"}, \code{"SONSCA"}, or \code{"DONSCA"}. Otherwise it is NULL.}
\item{tau}{The index of Goodman and Kruskal is given when the input parameter  \code{catype} is \code{"NSCA"}, \code{"SONSCA"}, or \code{"DONSCA"}. Otherwise it is  NULL.}
\item{inertiasum}{The total inertia of the analysis based on Pearson's chi-squared when catype is  \code{"CA"}, \code{"DOCA"} or \code{"SOCA"},
or based on the Goodman-Kruskal tau when \code{catype} is  \code{"NSCA"}, \code{"DONSCA"} or \code{"SONSCA"}  
(numerator of the Goodman-Kruskal tau index).}
\item{singvalue}{The singular values of the two-way contingency table.}
\item{inertias}{The inertia in absolute value and percentage, in the row space for each principal or polynomial axis.}
\item{inertias2}{The inertia in absolute value and percentage, in the column space for each principal or polynomial axis. 
When \code{catype} is \code{"CA"} or \code{"NSCA"} the associated inertia in the row and column spaces are the same for each principal axis.}
\item{t.inertia}{The total inertia of the two-way contingency table.}
\item{comps}{The polynomial components of inertia when the variables are ordered.}
\item{catype}{The type of correspondence analysis chosen by the analyst. By default, \code{catype = "ca"}.}
\item{mj}{The ordered scores of the column variable. When \code{mj = NULL}, the natural scores are used (\code{i = 1,...,cols}).}
\item{mi}{The ordered scores of the row variable. When \code{mi = NULL}, the natural scores are used (\code{i = 1,...,rows}).}
\item{pcc}{The weighted centered column profile matrix.}
\item{Jmass}{The weight matrix of the column variable.}
\item{Imass}{The weight matrix of the row variable.}
\item{Innprod}{The inner product, \code{Inner product}, of the biplot coordinates (for the two axes defined by \code{firstaxis = 1} and \code{lastaxis = 2})
}
\item{Z}{The generalised correlation matrix  when \code{catype = "SOCA"}, \code{catype = "DOCA"} , \code{catype = "SONSCA"}, \code{catype = "DONSCA"},  but
when \code{catype = "CA"}, or  \code{catype = "NSCA"}, it gives again the inner product matrix of biplot coordinates.
}
\item{M}{
The number of axes used for determining the structure of the elliptical confidence regions. 
By default,  \code{M = min(nrow(Xtable), ncol(Xtable)) - 1}, i.e. the rank of the data matrix.
}
%\item{ellcomp}{This parameter specifies whether the characteristics of the confidence ellipses (eccentricity, semi-axis, area, p-values) 
%are to be computed. By default, \code{ellcomp =  TRUE}.  
%}
\item{eccentricity}{
When  \code{ellcomp = TRUE}, the output gives  the eccentricity of the confidence ellipses.  
}
\item{row.summ}{
When  \code{ellcomp = TRUE}, the output gives for each row the summary results that contain the 
semi-major axis length of the ellipse, \code{HL Axis 1},  the semi-minor axis length for the ellipse,  \code{HL Axis 2}, 
the area of the ellipse, \code{Area} and the p-value, \code{P-value}. 
}
\item{col.summ}{
When  \code{ellcomp = TRUE}, the output gives for each column point, the summary results that contain the 
semi-major axis length of the ellipse, \code{HL Axis 1},  the semi-minor axis length for the ellipse,  \code{HL Axis 2}, 
the area of the ellipse, \code{Area} and the p-value, \code{P-value}. 
}
}

\references{
Beh EJ and Lombardo R 2014 Correspondence Analysis: Theory, Practice and New Strategies. Wiley.\cr
Lombardo R Beh EJ 2016 Variants of Simple Correspondence Analysis. The R Journal, 8 (2), 167--184.\cr
Lombardo R Beh EJ and Kroonenberg PM 2016 Modelling Trends in Ordered Correspondence Analysis Using Orthogonal
Polynomials. Psychometrika, 81(2), 325--349.
}
\author{
Rosaria Lombardo and Eric J Beh
}
\note{
This function recalls internally many other functions, depending on the setting of the input parameter \code{catype}, it recalls
one of the six functions which does a variant of correspondence analysis.
After performing a variant of correspondence analysis, it gives the output object necessary for printing and plotting the results. These two
important functions are \code{print.CAvariants} and \code{plot.CAvariants}. 
}
\examples{
data(asbestos)
CAvariants(asbestos, catype = "CA") 
CAvariants(asbestos, catype = "DOCA", mi = c(1:nrow(asbestos)), mj =c(4.5,14.5,24.5,34.5,44.5),  
firstaxis = 1, lastaxis = 2, M = min(nrow(asbestos), ncol(asbestos)) - 1) 
CAvariants(asbestos, catype = "DONSCA") 
data(shopdataM)
CAvariants(shopdataM, catype = "NSCA")
CAvariants(shopdataM, catype = "SONSCA")
CAvariants(shopdataM, catype = "SOCA")
}
\keyword{correspondence analysis variants}
\keyword{multivariate}% __ONLY ONE__ keyword per line
