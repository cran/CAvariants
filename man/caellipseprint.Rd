\name{caellipseprint}
\alias{caellipseprint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Numerical characteristics of elliptical confidence regions in variants of symmetrical correspondence analysis 
}
\description{
It produces the numerical characteristics of the elliptical confidence regions 
when symmetrical or ordered symmetrical correspondence analysis is performed. 
This function allows the analyst to superimpose confidence ellipses onto a graphical display when  \code{catype = "CA"}, \code{catype = "SOCA"} or \code{catype = "DOCA"}. 
It is called internally from the main plot function \code{print.CAvariants}. 
It uses the function \code{ellipse}.
}
\usage{
caellipseprint(Xtable, a1 = 1, a2 = 2, alpha = 0.05,  M = 2,  
Imass, Jmass, a, b, g, f, dmu, inertiapc)
}
\arguments{
\item{Xtable}{
The two-way contingency table.
}
\item{a1}{
The first axis number.
}
  \item{a2}{
The second axis number.
}
  \item{alpha}{
The confidence level of the elliptical regions. By default, \code{alpha = 0.05}.
}
  \item{M}{
The number of  axes used for constructing the ellypses. By default,  \code{M = 2}. Its maximum value is equal to the rank of the data matrix.}
   \item{Imass}{
The weight matrix of the row variable.
}
  \item{Jmass}{
The weight matrix of the column variable.
}
  \item{a}{
The row standard coordinates or, in case of the ordered variants of CA, the row standard polynomial coordinates.
}
  \item{b}{
The column standard coordinates or, in case of the ordered variants of CA, the column standard polynomial coordinates.
}
  \item{f}{
The row principal coordinates (scaled by a constant, by default \code{scaleplot = 1}).
}
  \item{g}{
The column principal coordinates (scaled by a constant, by default \code{scaleplot = 1}).
}
  \item{dmu}{
The squared singular values or inertia of each axis.
}
  \item{inertiapc}{
The percentage of explained inertia.}
}
\details{
The output values of this function.
}
\value{
\item{eccentricity}{Value of ellipse eccentricity. This is the distance between the ellipse center and either of its two foci, 
which can be thought of as a measure of how much the conic section deviates from being circular (when  it is equal to zero then the region becomes circular). }
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
Rosaria Lombardo and Eric J Beh
}
\keyword{ multivariate }% __ONLY ONE__ keyword per line
