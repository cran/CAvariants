\name{nscaellipseprint}
\alias{nscaellipseprint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Numerical characteristics of elliptical confidence regions in variants of non symmetrical correspondence analysis 
}
\description{
This function produces numerical characteristics of elliptical confidence regions 
when non symmetrical or ordered non symmetrical correspondence analysis is performed. 
It  superimposes the confidence ellipses on to the graphical displays, when 
 \code{catype = "NSCA"}, \code{catype = "SONSCA"} or \code{catype = "DONSCA"}.
It is called  from the main plot function \code{print.CAvariants} and uses the function \code{ellipse}.
}
\usage{
nscaellipseprint(Xtable, a1 = 1, a2 = 2, alpha = 0.05, M = 2, 
Imass, Jmass, a, b, f, g, dmu, tauden, inertiapc)
}
\arguments{
\item{Xtable}{
The two-way contingency table.
}
\item{a1}{
The dimension reflected along the horizontal axis.
}
  \item{a2}{
The dimension reflected along the vertical axis.
}
  \item{alpha}{
The confidence level of the elliptical regions. By default, \code{alpha = 0.05}.
}
   \item{M}{
The number of  axes used when drawing ellypse. By default,  \code{M = 2}.
Its maximum value is equal to the rank of the data matrix.}
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
  \item{f}{
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
The percentage of explained inertia.}
}
\details{
The output values of this function
}
\value{
  \item{eccentricity}{Value of ellipse eccentricity. This is the distance between the ellipse center and either of its two foci, 
which can be thought of as a measure of how much the conic section deviates from being circular (when  it is equal to zero then the region becomes circular).}
  \item{HL Axis 1}{Value of ellipse semi-axis 1 for each row and column points.}
  \item{HL Axis 2}{Value of ellipse semi-axis 2 for each row and column points.}
 \item{Area }{Ellipse area  for each row and column points.}
  \item{pvalcol }{P-value for each row and column points.}
}
\references{
Beh EJ and Lombardo R 2014 Correspondence Analysis: Theory, Practice and New Strategies. Wiley.\cr
Lombardo R Beh EJ 2016 Variants of Simple Correspondence Analysis. The R Journal, 8 (2), 167--184.
}
\author{
Rosaria Lombardo and Eric J Beh
}
\keyword{multivariate}% __ONLY ONE__ keyword per line
