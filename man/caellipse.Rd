\name{caellipse}
\alias{caellipse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Algebraic elliptical confidence regions 
}
\description{
It produces elliptical confidence regions when symmetrical or ordered symmetrical correspondence analysis is performed. \cr
This function allows the analyst to superimpose confidence ellipses into a graphical display when the input parameter \code{catype}
of the main code \code{CAvariants} is equal to \code{CA, SOCA} or \code{DOCA}. \cr
It is called internally  from the main plot function \code{plot.CAvariants}. \cr
It uses the function \code{ellipse}.
}
\usage{
caellipse(N, a1 = 1, a2 = 2, alpha = 0.05, cols = c(2, 4), 
M = 2, cex = 0.8, cex.lab = 0.5,  
mar = c(5, 4, 4, 2) + 0.1, prop = 0.8, Imass, 
Jmass, a, b, g, fr, dmu, inertiapc,  plottype="biplot", 
biptype = "row",pos=2,  arrow = T, length = 0,graphy=T)
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
The graphical parameter for setting the colours of points in the graphical displays.
}
  \item{M}{
The number of  axes used when drawing ellypse. By default, it is equal to \code{M= 2}. Its maximum value is equal to the rank of the data matrix.}
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
  \item{inertiapc}{
The percentage of explained inertia.
}
  \item{pos}{
The parameter for specifying the position of  point symbols in the graphical displays. \cr
By default, \code{pos=2}.}
  \item{plottype}{
The parameter for specifying the type of plot. By default, \code{plottype="classic"}. 
}
 \item{biptype}{
The parameter for specifying the type of biplot. By default, \code{biptype="row"}. 
}
\item{arrow}{
The parameter used for diplaying arrows in biplot in correspondence of points in principal coordinates. By default, \code{arrow=TRUE}. 
}
\item{length}{The parameter used for setting the length of the arrow end in biplot. By default, \code{length=0}. }
\item{graphy}{The flag parameter used for  diplaying ellipses when using \code{plot.CAvariants} and not displaying when using \code{print.CAvariants} function.
 It is set as \code{graphy=T} in \code{plot.CAvariants}  and as \code{graphy=F}  in \code{print.CAvariants}.}
}
\details{
The output values of this function
}
\value{
\item{eccentricity}{Value of ellipse eccentricity.}
\item{HL Axis 1}{Value of ellipse semi-axis 1 for each row and column points.}
\item{HL Axis 2}{Value of ellipse semi-axis 2 for each row and column points.}
\item{Area}{Ellipse area  for each row and column points.}
\item{pvalcol}{P-value for each row and column points.}
}
\references{
Beh EJ and Lombardo R 2014 Correspondence Analysis: Theory, Practice and New Strategies. John Wiley & Sons.
}
\author{
Rosaria Lombardo and Eric J Beh
}
\note{
This function is called from the main plot function \code{plot.CAvariants} and is executed when 
\code{catype="CA"}, \code{catype="SOCA"} or \code{"DOCA"} only if \code{ell=T}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\keyword{ multivariate }% __ONLY ONE__ keyword per line
