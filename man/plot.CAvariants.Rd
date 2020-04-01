\name{plot.CAvariants}
\alias{plot.CAvariants}
\title{
Main plot function 
}
\description{
This function produces the graphical display for the selected variant of correspondence analysis. 
When \code{catype = "CA"} \code{catype = "NSCA"} and \code{plottype = "classic"}, the function produces a plot 
of the principal coordinates for the row and column categories.\cr
When \code{plottype = "biplot"}, it produces a biplot graphical display, or a polynomial biplot
in case of ordered variables.
For an ordered analysis only the polynomial biplots are constructed. In particular, 
for the singly ordered variants only the row isometric polynomial biplot is appropriate.
When the parameter \code{catype} defines an ordered variant of CA, the input parameter  
\code{plottype} should be equal to \code{plottype = "biplot"}.  If \code{biptype = "row"}, 
it will produce a row isometric polynomial biplot.
}
\usage{

\S3method{plot}{CAvariants}(x, firstaxis = 1, lastaxis = 2, thirdaxis = 3, cex = 0.8, 
cex.lab = 0.8, plottype = "biplot", biptype = "row",  
scaleplot = 1, posleg = "topleft", pos = 2, ell = FALSE,  
alpha = 0.05, plot3d = FALSE, size1 = 1.5, size2 = 3, invproj = TRUE,\dots)
}
\arguments{
  \item{x}{
The name of the output object used with the main function \code{CAvariants}.
}
  \item{firstaxis}{
The horizontal polynomial, or principal, axis. By default, \code{firstaxis = 1}.
}
  \item{lastaxis}{
The vertical polynomial, or principal, axis. By default, \code{lastaxis = 2}.
}
 \item{thirdaxis}{
The third polynomial, or principal, axis in tridimensional plot. By default, \code{thirdaxis = 3}.
}
  \item{cex}{
The parameter for setting the size of the character labels for the points in a graphical display.  By default,  \code{cex = 0.8}.
}
 \item{cex.lab}{
The parameter for setting the size of the character labels of axes in graphical displays.  By default,  \code{cex.lab = 0.8}.
}
  \item{plottype}{
The type of graphical display required (either a correspondence plot or a biplot). 
The type of graphical display to be constructed. By default, \code{plottype = "biplot"}; 
the alternative is \code{plottype = "classic"}. 
}
\item{biptype}{
The parameter for specifying the type of biplot. By default, \code{biptype = "row"}. 
One may specify a row-isometric biplot  (\code{biptype = "row"}) or a column-isometric biplot (\code{biptype = "column"}).  
This feature is available for the nominal symmetrical and non-symmetrical correspondence analyses. 
By default, a row-isometric biplot, \code{biptype = "row"}, is produced.
}
\item{scaleplot}{
The parameter for scaling the classic plot and biplot coordinates. See  Gower et al. (2011), section 2.3.1, or 
page 135 of Beh and Lombardo (2014). By default, \code{scaleplot = 1}.
}
\item{posleg}{
The position of the legend when portraying trends of the categories 
for ordered variants of correspondence analysis. 
By default, \code{posleg = "topleft"}.}
 \item{pos}{
The parameter that specifies the position of label of each point in the graphical display. By default, \code{pos = 2}.
}
\item{ell}{The logical parameter which specifies whether algebraic confidence ellipses are to be included in the plot or not.  
Setting the input parameter to \code{ell = TRUE} will assess the statistical significance of each category to 
the association between the variables. By default,  \code{ell = FALSE}.}
%\item{M}{
%The number of axes  considered in determining the structure of the elliptical confidence regions.  
%By default,  \code{M = min(nrow(Xtable), ncol(Xtable)) - 1}, i.e. the rank of the data matrix.
%}
 \item{alpha}{
The confidence level of the elliptical regions. By default, \code{alpha = 0.05}.
}
\item{plot3d}{The logical parameter specifies whether a 3D plot is to be included 
in the output or not. By default, \code{plot3d = FALSE}.
}
\item{size1}{
The size of the plotted symbol. By default, \code{size = 1.5}.
}
\item{size2}{
The size of the plotted text. By default, \code{size = 3}.
}
\item{invproj}{The logical parameter specifies whether to portray standard coordinates as vectors and principal coordinates as points or vice-versa.
By default, \code{invproj = TRUE}.}
\item{\dots}{Further arguments passed to, or from, other functions.}
}
\details{
It produces either a classical or biplot graphical display. Further, when \code{catype = "DOCA"},
\code{catype = "SOCA"}, \code{catype = "DONSCA"} or \code{catype = "SONSCA"}, 
the trends of the row and column variables (after the reconstruction of column profiles by the polynomials) is portrayed.\cr
For classical biplot displays, it superimposes the algebraic confidence ellipses. It uses the secondary plot function \code{caellipse} (or
\code{nscaellipse}) for the symmetrical (or non symmetrical) CA variants.}
\references{
Beh EJ and Lombardo R 2014 Correspondence Analysis: Theory, Practice and New Strategies. Wiley.\cr
Gower J, Lubbe S, and le Roux, N 2011 Understanding Biplots. Wiley.\cr
Lombardo R Beh EJ 2016 Variants of Simple Correspondence Analysis. The R Journal, 8 (2), 167--184.\cr
Lombardo R Beh EJ and Kroonenberg PM 2016 Modelling Trends in Ordered Correspondence Analysis Using Orthogonal
Polynomials. Psychometrika, 81(2), 325--349.
}
\author{
Rosaria Lombardo and Eric J Beh
}
\note{
For the classical plots, row and column principal coordinates are plotted.
For biplots, one set of coordinates is the standard coordinates and the other 
is the principal coordinates. 
When an ordered variant of correspondence analysis is performed,
the biplot is constructed where one set of coordinates consists of 
the standard polynomial coordinates and the other one is the principal polynomial coordinates. 
}
\examples{
data(asbestos)
resasbestosCA<-CAvariants(asbestos, catype = "CA", M=2) 
plot(resasbestosCA, plottype = "classic", plot3d = TRUE)
plot(resasbestosCA, plottype = "classic",  ell = TRUE)
plot(resasbestosCA, plottype = "biplot", biptype = "column", scaleplot=1.5, invproj=FALSE)
resasbestosDOCA<-CAvariants(asbestos, catype = "DOCA") 
plot(resasbestosDOCA, plottype = "biplot", biptype = "column")
resasbestosNSCA<-CAvariants(asbestos, catype = "NSCA") 
plot(resasbestosNSCA, plottype = "biplot", biptype = "column", plot3d = TRUE)
}
\keyword{multivariate}% __ONLY ONE__ keyword per line
