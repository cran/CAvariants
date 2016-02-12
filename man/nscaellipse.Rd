\name{nscaellipse}
\alias{nscaellipse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Algebraic elliptical confidence regions in non symmetrical correspondence analysis 
}
\description{
This function produces elliptical confidence regions when non symmetrical or ordered non symmetrical correspondence analysis is performed. 
It allows to superimpose confidence ellipses in graphical displays, when the input parameter \code{catype} 
of the main code \code{CAvariants} is equal to \code{"NSCA", "SONSCA"} or \code{"DONSCA"}.\cr
It is called  from the main plot function \code{plot.CAvariants} and uses the function \code{ellipse}.
}
\usage{
nscaellipse(N, a1 = 1, a2 = 2, alpha = 0.05, cols = c(2, 4), 
M = 2, cex = 0.8, cex.lab = 0.5, 
mar = c(5, 4, 4, 2) + 0.1, prop = 0.8,  Imass, Jmass,
a, b, fr, g, dmu,tauden, inertiapc, plottype = "biplot", 
biptype = "row", pos = 2, arrow = T, length = 0,graphy=T)
}
\arguments{
\item{N}{
The two-way contingency table.
}
\item{a1}{
The dimension reflected along the horizontal axis.
}
  \item{a2}{
The dimension reflected along the vertical axis.
}
  \item{alpha}{
The confidence level of the elliptical regions.
}
  \item{cols}{
The graphical parameter for setting the colours of points in graphical displays.
}
  \item{M}{
The number of  axes used when drawing ellypse. By default, it is equal to \code{M= 2}.
Its maximum value is equal to the rank of the data matrix.}
  \item{cex}{
The parameter for setting the size of character labels of points in graphical displays.
}
  \item{cex.lab}{
The parameter for setting the size of character labels of axes in graphical displays.
}
  \item{mar}{
The parameter for setting the size of the plotting area. 
}
  \item{prop}{
The scaling parameter for specifying the limits of the plotting area. By default it is equal to 1. 
}
  \item{Imass}{
The weight matrix of the row variable.
}
  \item{Jmass}{
The weight matrix of the column variable.
}
  \item{a}{
The row principal or polynomial axes.
}
  \item{b}{
The column principal or polynomial axes.
}
  \item{fr}{
The row coordinates.
}
  \item{g}{
The column coordinates.
}
  \item{dmu}{
The squared singular values or inertia of each axis.
}
  \item{tauden}{The denominator of tau index.
}
\item{inertiapc}{
The percentage of explained inertia.
}
\item{plottype}{
The kind of graphical display. By default, \code{plottype="biplot"}.
}
\item{biptype}{
The parameter for specifying the type of biplot. By default, \code{biptype="row"}. 
}
\item{pos}{
The parameter for specifying the position of the point symbols in the graphical displays. 
By default, \code{pos=2}.
}
\item{arrow}{
The parameter for diplaying arrows in biplot in correspondence of points in principal coordinates. By default, \code{arrow=TRUE}. 
}
\item{length}{
The parameter for setting the length of the arrow end in biplot. By default, \code{length=0.01}. 
}
\item{graphy}{The parameter for diplaying ellipses in plots or  biplots. By default, \code{graphy=TRUE}.}
}
\details{
The output values of this function
}
\value{
  \item{eccentricity}{Value of ellipse eccentricity.}
  \item{HL Axis 1}{Value of ellipse semi-axis 1 for each row and column points.}
  \item{HL Axis 2}{Value of ellipse semi-axis 2 for each row and column points.}
 \item{Area }{Ellipse area  for each row and column points.}
  \item{pvalcol }{P-value for each row and column points.}

}
\references{
Beh EJ and Lombardo R 2014 Correspondence Analysis: Theory, Practice and New Strategies. John Wiley & Sons.
}
\author{
Rosaria Lombardo and Eric J Beh
}
\note{
This function is called from the main plot function \code{plot.CAvariants} and 
is executed when \code{catype="NSCA"}, \code{catype="SONSCA"} or 
\code{catype="DONSCA"}.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%\examples{
%  }
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{multivariate}
\keyword{ multivariate }% __ONLY ONE__ keyword per line
