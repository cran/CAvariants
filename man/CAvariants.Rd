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
CAvariants(
Xtable, mj = NULL, mi = NULL, firstaxis = 1, lastaxis = 2,
catype = "CA", ellcomp = TRUE, Mell = min(nrow(Xtable), ncol(Xtable)) - 1, alpha = 0.05) 
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
The horizontal polynomial or principal axis.  By default \code{firstaxis = 1}.
}
\item{lastaxis}{
The vertical polynomial or principal axis.  By default \code{lastaxis = 2}.
}
\item{catype}{
The input parameter for specifying what variant of correspondence analysis is considered. By default,  \code{catype = "CA"}. 
Other possible values are:  \code{catype = "SOCA"}, \code{catype = "DOCA"}, \code{catype = "NSCA"}, \code{catype = "SONSCA"}, \code{catype = "DONSCA"}.
}
\item{ellcomp}{This input parameter ensures that the characteristics of the algebraic confidence ellipses are computed and stored. 
When \code{ellcomp = TRUE} (which is by default), the output includes the characteristics of the ellipses. 
The eccentricity of the confidence ellipses is summarised by the quantity \code{eccentricity}, this is the distance between the center and either 
of its two foci,
which can be thought of as a measure of how much the conic section deviates from being circular
 
(when  it is equal to zero then the region becomes circular).
 The semi-major axis length of the ellipse for each row and column 
point is given by \code{HL Axis 1} while \code{HL Axis 2} gives the semi-minor axis length of the points along the second axis. 
The area of the ellipse for each row and column category is given by \code{Area} while the p-value of each category  is defined by \code{P-value}.}
\item{Mell}{
The number of axes \code{Mell} considered in determining the structure of the elliptical confidence regions. 
By default,  \code{Mell = min(nrow(Xtable), ncol(Xtable)) - 1}, i.e. the rank of the data matrix.
}
\item{alpha}{
The confidence level of the elliptical regions. By default, \code{alpha = 0.05}.
}
}
\details{
This function belongs to the object class called \code{cacorporateplus}  
}
\value{
Description of the output returned
\item{Xtable }{The two-way contingency table.}
\item{rows}{The row number of the two-way contingency table.}
\item{cols}{The column number of the two-way contingency table.}
\item{r}{The rank of the two-way contingency table.}
\item{rowlabels}{The label of the row variable.}
\item{collabels}{The label of the column variable.}
\item{Rprinccoord}{The row principal  coordinates.  When the input parameter \code{catype} is \cr
\code{"DOCA"}, \code{"SOCA"}, \code{"SONSCA"} or \code{"DONSCA"}, they are row principal polynomial coordinates.}
\item{Cprinccoord}{The column principal  coordinates. When the input parameter \code{catype} is  \cr
\code{"DOCA"} or \code{"DONSCA"}, they are column principal polynomial coordinates.}
\item{Rstdcoord }{The row standard  coordinates.  When the input parameter \code{catype} is \cr
\code{"DOCA"} or \code{"DONSCA"},  they are row standard polynomial coordinates.}
\item{Cstdcoord }{The column standard  coordinates. When the input parameter  \code{catype} is  \cr
\code{"DOCA"}, \code{"SOCA"}, \code{"SONSCA"} or \code{"DONSCA"},  they are column standard polynomial coordinates.}
\item{tauden}{The tau denominator is given when the input parameter  \code{catype}  is \code{"NSCA"}, \code{"SONSCA"}, or \code{"DONSCA"}, otherwise it is a null value.}
\item{tau}{The tau index is given when the input parameter  \code{catype} is \code{"NSCA"}, \code{"SONSCA"}, or \code{"DONSCA"}, otherwise it is a null value.}
\item{inertiasum}{The total inertia of the classical correspondence analysis when catype is  \code{"CA"}, \code{"DOCA"} or \code{"SOCA"} (the Pearson's index), 
or the inertia of non symmetrical correspondence analysis when \code{catype} is  \code{"NSCA"}, \code{"DONSCA"} or \code{"SONSCA"}  
(numerator of the Goodman-Kruskal tau index).}
\item{inertias}{The associated inertia in absolute value and percentage, in the row space for each principal or polynomial axis.}
\item{inertias2}{The associated inertia in absolute value and percentage, in the column space for each principal or polynomial axis. 
When \code{catype} is \code{"CA"} or \code{"NSCA"} the associated inertia in the row and column spaces are the same for each principal axis.}
\item{comps}{The polynomial components of inertia when the variables are ordered. The inertia of row and/or column space is partitioned in terms 
of polynomial components in ordered CA variants.}
\item{catype}{The kind of correspondence analysis chosen.}
\item{mj}{The ordered scores of a column variable. When \code{mj = NULL}, the natural ordered numbers are shown.}
\item{mi}{The ordered scores of a row variable. When \code{mi = NULL}, the natural ordered numbers are shown.}
\item{pcc}{The weighted centered column profile matrix.}
\item{Jmass}{The weight matrix of the column variable.}
\item{Imass}{The weight matrix of the row variable.}
\item{Trend}{The inner product, \code{Inner product}, of the biplot coordinates (concerning the first two axes when \code{firstaxis=1} and \code{lastaxis=2})
}
\item{Z}{The generalized correlation matrix  when catype is  \code{"SOCA"}, \code{"DOCA"} , \code{"SONSCA"}, \code{"DONSCA"},  but
when catype is  \code{"CA"}, \code{"NSCA"}, it gives again the inner product matrix of biplot coordinates.
}
\item{ellcomp}{The flag parameter, \code{ellcomp}, specifies that the characteristics of the confidence ellipses (eccentricity, semi-axis, area, p-values) 
are computed. \cr
By default, \code{ellcomp =  TRUE}.  
}
\item{risell}{
When the input parameter, \code{ellcomp}, is set to \code{ellcomp = TRUE}, 
the output includes the  characteristics \code{risell} of the confidence ellipses, the eccentricity of the confidence ellipses, \code{risell$eccentricity},  
for each row and column point, the summary results,  \code{risell$row.summ} and \code{risell$col.summ}, contain the 
semi-major axis length of the ellipse, \code{HL Axis 1},  the semi-minor axis length for the ellipse,  \code{HL Axis 2}, 
the area of the ellipse, \code{Area} and the p-value, \code{P-value}. 
}
\item{Mell}{The number of axes \code{Mell} considered in determining the structure of the elliptical confidence regions.  
By default,  \code{Mell = min(nrow(Xtable), ncol(Xtable)) - 1}, i.e. the rank of the data matrix.
}
}

\references{
Beh EJ and Lombardo R 2014 Correspondence Analysis: Theory, Practice and New Strategies. John Wiley & Sons.
Lombardo R Beh EJ Variants of Simple Correspondence Analysis. The R Journal (accepted).
}
\author{
Rosaria Lombardo and Eric J Beh
}
\note{
This function recalls internally many other functions, depending on the setting of the input parameter \code{catype}, it recalls
one of the six functions which does a variant of correspondence analysis.
After performing a variant of correspondence analysis, it gives the output object necessary for printing and plotting the results. These two
important functions are \code{print.CAvariants} and \code{plot.CAvariants}. 
This function belongs to the class \code{cacorporateplus}. 
}
\examples{
data(asbestos)
CAvariants(asbestos, catype = "CA", firstaxis = 1, lastaxis = 2) 
CAvariants(asbestos, catype = "DOCA", firstaxis = 1,lastaxis = 2) 
CAvariants(asbestos, catype = "DONSCA",firstaxis=1, lastaxis = 2, ellcomp = FALSE) 
data(shopdataM)
CAvariants(shopdataM, catype = "NSCA", firstaxis = 1, lastaxis = 2)
CAvariants(shopdataM, catype = "SONSCA", firstaxis = 1, lastaxis = 2)
CAvariants(shopdataM, catype = "SOCA", firstaxis = 1, lastaxis = 2)
}
\keyword{correspondence analysis variants}
\keyword{multivariate}% __ONLY ONE__ keyword per line
