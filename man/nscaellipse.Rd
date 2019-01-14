\name{nscaellipse}
\alias{nscaellipse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Algebraic elliptical confidence regions for non-symmetrical variants of correspondence analysis 
}
\description{
This function produces elliptical confidence regions when non-symmetrical or ordered non-symmetrical correspondence analysis is performed. 
It  superimposes the confidence ellipses onto the graphical displays, when \code{catype = "NSCA"}, 
\code{catype = "SONSCA"} or \code{catype = "DONSCA"}.
It is called  from the main plot function \code{plot.CAvariants} and uses the function \code{ellipse}.
}
\usage{
nscaellipse(Xtable, a1 = 1, a2 = 2, alpha = 0.05, cols = c(2, 4), 
M = 2, cex = 0.8, cex.lab = 0.8, mar = c(5, 4, 4, 2) + 0.1, prop = 0.8,  Imass, Jmass,
a, b, f, g, dmu, tauden, inertiapc, plottype = "biplot", 
biptype = "row", pos = 2, arrow = TRUE, length = 0,  ell = TRUE)
}
\arguments{
\item{Xtable}{
The two-way contingency table.
}
\item{a1}{
The axis number of the horizontal axis.
}
  \item{a2}{
The axis number of the vertical axis.
}
  \item{alpha}{
The confidence level of the elliptical regions. By default, \code{alpha = 0.05}.
}
  \item{cols}{
The graphical parameter for setting the colours of the points in the graphical displays.
}
  \item{M}{
The number of axes used for constructing the confidence ellypses. By default,  \code{M = 2}.
Its maximum value is equal to the rank of the data matrix.}
  \item{cex}{
The parameter for setting the size of the character labels for the points in a graphical display.  By default,  \code{cex = 0.8}.
}
 \item{cex.lab}{
The parameter for setting the size of the character labels of axes in graphical displays.  By default,  \code{cex.lab = 0.8}.
}
 \item{mar}{
The parameter for setting the size of the plotting area. By default, \code{mar = c(5, 4, 4, 2) + 0.1}.
}
  \item{prop}{
The scaling parameter for specifying the limits of the plotting area. By default,  \code{prop = 1}. 
}
  \item{Imass}{
The weight matrix of the row variable.
}
  \item{Jmass}{
The weight matrix of the column variable.
}
  \item{a}{
The row standard coordinates or, in case of the ordered variants of NSCA, the row standard polynomial coordinates.
}
  \item{b}{
The column standard coordinates or, in case of the ordered variants of NSCA, the column standard polynomial coordinates.
}
  \item{f}{
The row principal coordinates (scaled by a constant, by default \code{scaleplot = 1}).
}
  \item{g}{
The column principal coordinates (scaled by a constant, by default \code{scaleplot = 1}).
}
  \item{dmu}{
The squared singular values, or principal inertia, of each axis.
}
  \item{tauden}{The denominator of the tau index.
}
\item{inertiapc}{
The percentage of explained inertia along each of the axes.
}
\item{plottype}{
The type of graphical display to be constructed. By default, \code{plottype = "biplot"}; the alternative is \code{plottype = "classic"}. 
}
\item{biptype}{
The parameter for specifying the type of biplot. By default, \code{biptype = "row"}. 
}
\item{pos}{
The parameter that specifies the position of label of each point in the graphical display. By default, \code{pos = 2}.
}
\item{arrow}{
The parameter for diplaying the arrows in a biplot. By default, \code{arrow = TRUE}. 
}
\item{length}{
The parameter for setting the length of the arrow end in a biplot. By default, \code{length = 0}. 
}
\item{ell}{The logical parameter used for displaying the confidence ellipses. By default, \code{ell = TRUE}.}
}
\details{
The output values of this function
}
\value{
  \item{eccentricity}{The eccentricity of the ellipses. This is the distance between the centre of the ellipse and its two foci, 
which can be thought of as a measure of how much the conic section deviates from being circular (when the region is perfectly circular, eccentricity is zero). }
\item{HL Axis 1}{Value of the semi-major axis length for each row and column point.}
\item{HL Axis 2}{Value of the semi-minor axis length for each row and column point.}
\item{Area}{Area of the ellipse for each row and column point.}
\item{pvalcol}{Approximate p-value for each of the row and column points.}
}
\references{
Beh EJ and Lombardo R 2014 Correspondence Analysis: Theory, Practice and New Strategies. Wiley.\cr
Lombardo R Beh EJ 2016 Variants of Simple Correspondence Analysis. The R Journal, 8 (2), 167--184.
}
\author{
Rosaria Lombardo and Eric J Beh
}
\note{
This function is called from the main plot function \code{plot.CAvariants} and 
is executed when \code{catype = "NSCA"}, \code{catype = "SONSCA"} or 
\code{catype = "DONSCA"},  only if \code{ell = TRUE}.
}
\keyword{ multivariate }% __ONLY ONE__ keyword per line
