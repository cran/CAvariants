\name{plot.CAvariants}
\alias{plot.CAvariants}
\title{
Main plot function 
}
\description{
This function allows the analyst to produce the suitable graphical displays with respect to the six variants of correspondence analysis. 
In particular when \code{plottype = "classic"}, it produces classical graphical displays for \code{catype = "CA"} and \code{catype = "NSCA"}, 
where the row and column variables are graphically depicted in principal coordinates.\cr
When we set \code{plottype = "biplot"}, it produces biplot graphical displays, or polynomial biplots in case of ordered analysis.
Notice that for ordered analysis only polynomial biplots are suitable.  In particular, for the singly ordered variants only row isometric polynomial biplots make sense,
as we assume that the ordered variable is the column variable (the column coordinates are standard polynomial coordinates 
and the row coordinates are principal polynomial coordinates).
When the input parameter \code{catype}  concerns an ordered variant of CA, then the input parameter  
\code{plottype} should be equal to \code{plottype = "biplot"},  if \code{biptype = "row"}, 
it will give back a row isometric polynomial biplot.
 }
\usage{
\S3method{plot}{CAvariants}(x, firstaxis = 1, lastaxis = 2, cex = 0.8, 
cex.lab = 0.8, prop = 1,  plottype = "biplot", biptype = "row",  
scaleplot = 1, posleg = "topleft", pos = 2, ell = FALSE, Mell = x$Mell, 
alpha = 0.05, plot3d = TRUE, size=5,adj=c(0,0.5),\dots)
}
\arguments{
  \item{x}{
The name of the output object, for example say \code{res}, used with the main function \code{CAvariants}.
}
  \item{firstaxis}{
The horizontal polynomial or principal axis, \code{firstaxis}. By default, \code{firstaxis = 1}.
}
  \item{lastaxis}{
The vertical polynomial or principal axis, \code{lastaxis}. By default, \code{lastaxis = 2}.
}
  \item{cex}{
The size of characters, \code{cex}, displayed on the correspondence plot or biplot. By default, \code{cex = 0.8}.}
 \item{cex.lab}{
The parameter \code{cex.lab} that specifies the size of character labels of axes in graphical displays. By default, \code{cex.lab = 0.8}.
}
  \item{prop}{
The scaling parameter for specifying the limits of the plotting area. By default, \code{prop = 1}. 
}
 \item{plottype}{
The type of graphical display required (either a classical correspondence plot or a biplot). 
The user can look at a classical correspondence plot by defining the input
parameter  \code{plottype = "classic"}.
When \code{plottype = "biplot"}, it produces biplot graphical displays, or polynomial biplots in case of an ordered analysis.

Note that for ordered analysis only polynomial biplots are suitable.  In particular  for the singly ordered variants, 
only row isometric polynomial biplots make sense, as we assume that the ordered variable is the column variable 
(the column coordinates are standard polynomial coordinates and the row coordinates are principal polynomial coordinates). 
 By default, \code{plottype = "biplot"}.
}
\item{biptype}{
For a biplot, one may specify that it be a row-isometric biplot  (\code{biptype = "row"}) or a column-isometric biplot (\code{biptype = "column"}).  
This feature is available for the nominal symmetrical and the non symmetrical correspondence analyses. 
By default, a row-isometric biplot, \code{biptype = "row"}, is produced.
}
\item{scaleplot}{
The parameter for scaling the biplot coordinates, \code{scaleplot}, originally proposed in Section 2.3.1 of Gower et al. (2011) and 
described on page 135 of  Beh and Lombardo (2014). By default, \code{scaleplot = 1}.
}
\item{posleg}{
The parameter \code{posleg} for specifying the position of the legend when portraying trends of ordered categories 
in ordered variants of correspondence analysis. 
By default, \code{posleg = "topleft"}. }
 \item{pos}{
The parameter for specifying the position of  point symbols in the graphical displays.
By default, \code{pos = 2}.}
\item{ell}{The logical parameter, \code{ell} which specifies whether algebraic confidence ellipses are to be included in the plot or not.  
Setting the input parameter to \code{ell = TRUE} will allow the user to assess the statistical significance of each category to 
the association between the variables. The ellipses will be included when the plot is constructed using principal coordinates 
(being either row and column principal coordinates or row and column principal polynomial coordinates). 
By default, this input parameter is set to \code{ell = FALSE}.  See also the input parameter \code{ellcomp}  
of the function \code{CAvariants} 
for a  description of the numeric characteristics of the confidence ellipses (eccentricity, area, etc.), as well as the
input parameter \code{ellprint} of the function \code{print.CAvariants} for getting a print of these parameters.
}
\item{Mell}{
The number of axes \code{Mell} considered when portraying  the elliptical confidence regions.  \cr
By default, it is equal to \code{Mell = min(nrow(Xtable), ncol(Xtable)) - 1}, i.e. the rank of the data matrix. 
This parameter is identical to the input parameter \code{Mell} of the function \code{CAvariants}.
}
 \item{alpha}{
The confidence level of the elliptical regions. By default, \code{alpha = 0.05}.}
\item{plot3d}{The logical parameter, \code{plot3d} which specifies whether 3D plot are to be included 
in the output or not. By default, \code{plot3d = FALSE}.
}
\item{size}{
The size of the plotting area.
}
\item{adj}{one or two values in [0, 1] which specify the x (and optionally y) adjustment of the labels.}
\item{\dots}{Further arguments passed to or from other methods.}
}
\details{
It produces classical and biplot graphical displays. Further when \code{catype} is equal to \cr
\code{"DOCA", "SOCA", "DONSCA"} or \code{"SONSCA"}, 
the trend of row and column variables after the reconstruction of column profiles by polynomials is portrayed.\cr
For classical biplot displays, it superimposes on it algebraic ellipses of confidence. It uses the secondary plot functions \code{caellipse} or
\code{nscaellipse}, depending on the input parameter \code{catype}. }
\references{
Beh EJ and Lombardo R 2014 Correspondence Analysis: Theory, Practice and New Strategies. John Wiley & Sons.\cr
Gower J, Lubbe S, and le Roux, N 2011 Understanding Biplots. John Wiley & Sons.\cr
Lombardo R Beh EJ 2016 Variants of Simple Correspondence Analysis. The R Journal, 8 (2), 167--184.\cr
Lombardo R Beh EJ and Kroonenberg PM 2016 Modelling Trends in Ordered Correspondence Analysis Using Orthogonal
Polynomials. Psychometrika, 81(2), 325--349.
}
\author{
Rosaria Lombardo and Eric J Beh
}
\note{
For classical graphical displays, both sets of coordinates are defined using principal coordinates (see Greenacre's terminology).
In biplot graphical displays, one set of coordinates is standard and the other is principal. 
When the analysis is ordered, it makes sense only biplot. One set of coordinates consists of standard polynomial coordinates and 
the other one is of principal polynomial coordinates.
}
\examples{
data(asbestos)
risasbestos<-CAvariants(asbestos, catype = "CA", firstaxis=1, lastaxis=2) 
plot(risasbestos, plottype = "classic", plot3d = TRUE)
plot(risasbestos, plottype = "biplot", biptype = "row", ell = TRUE)
}
\keyword{multivariate}% __ONLY ONE__ keyword per line
