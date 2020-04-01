\name{vcaellipse}
\alias{vcaellipse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Algebraic elliptical confidence regions for symmetrical variants of correspondence analysis 
}
\description{
It produces elliptical confidence regions when symmetrical or ordered symmetrical correspondence analysis is performed. 
This function allows the analyst to superimpose confidence ellipses onto a graphical display when the input parameter \code{catype}
of the main function \code{CAvariants} is set to \code{"CA", "SOCA"} or \code{"DOCA"}. 
It is called internally from the main plot function \code{plot.CAvariants}. 
It uses the function \code{ellipse}.
}
\usage{
vcaellipse(t.inertia, inertias, inertiapc, cord1, cord2, a, b, firstaxis=1,
lastaxis = 2, n, M = 2, Imass, Jmass)
}
\arguments{
\item{t.inertia}{
The total inertia of the two-way contingency table (Pearson's chi-squared or Goodman and Kruskal's index depends on the CA variant).
}
\item{inertias}{
The explained inertia of each dimension.
}
  \item{inertiapc}{
The percentage of explained inertia for each dimension.
}
  \item{cord1}{
The row principal coordinates. 
}
  \item{cord2}{
The column principal coordinates.
}
  \item{a}{
The row standard coordinates or, in case of the ordered variants of CA, the row standard polynomial coordinates.
}
  \item{b}{
The column standard coordinates or, in case of the ordered variants of CA, the column standard polynomial coordinates.
}
  \item{firstaxis}{
The horizontal polynomial, or principal, axis. By default, \code{firstaxis = 1}.
}
  \item{lastaxis}{
The vertical polynomial, or principal, axis. By default, \code{lastaxis = 2}.
}
  \item{n}{
The total number of observations. 
}
\item{M}{
The number of axes  considered in determining the structure of the elliptical confidence regions.  
}
   \item{Imass}{
The weight matrix of the row variable.
}
  \item{Jmass}{
The weight matrix of the column variable.
}
}
\details{
The output values of this function.
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
Beh EJ 2010 Elliptical confidence regions for simple correspondence analysis. J. Stat. Plan.
Inference 140, 2582--2588.\cr
Beh EJ and Lombardo R 2014 Correspondence Analysis: Theory, Practice and New Strategies. Wiley.\cr
Beh EJ Lombardo R 2015 Confidence regions and Approximate P-values for classical and non-symmetric correspondence analysis. 
Journal of Communications and Statistics, Theory and Methods. 44: 95--114.\cr
Lombardo R Beh EJ 2016 Variants of Simple Correspondence Analysis. The R Journal, 8 (2), 167--184.
}
\author{
Rosaria Lombardo and Eric J Beh
}
\note{
This function is called from the main plot function \code{plot.CAvariants} and is executed when in the main plot function the parameter
\code{ell = TRUE}.
}
\keyword{ multivariate }% __ONLY ONE__ keyword per line
