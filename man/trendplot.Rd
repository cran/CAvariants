\name{trendplot}
\alias{trendplot}
\title{
Plot of reconstructed centred row and column profiles 
}
\description{
This function portrays the trend of the centred row and column of the column profile matrix 
that are reconstructed using othogonal polynomials and/or principal axes.}
\usage{
trendplot(mf, mg, cex = 1, cex.lab = 0.8, prop = 0.5, posleg = "topleft", 
xlab = "First Principal Axis", ylab = "Second Principal Axis")
}
\arguments{
  \item{mf}{
The number of row or column categories. By defaults, natural numbers. 
}
  \item{mg}{
The inner product matrix.
}
  \item{cex}{
The parameter for setting the size of the character labels for the points in a graphical display.  By default,  \code{cex = 0.8}.
}
  \item{cex.lab}{
The parameter for setting the size of character labels of axes in graphical displays. By default, \code{cex.lab = 0.8}.
}
  \item{prop}{
The scaling parameter for specifying the limits of the plotting area. By default, \code{prop = 0.5}. 
}
  \item{posleg}{The position of the legend when portraying trends of the categories 
for ordered variants of correspondence analysis. 
By default, \code{posleg = "topleft"}.}
  \item{xlab}{
The parameter for setting the character label along the horizontal axis of the graphical display.
}
  \item{ylab}{
The parameter for setting the  character label along the vertical axis of the graphical display.
}
}
\references{
Beh EJ and Lombardo R 2014 Correspondence Analysis: Theory, Practice and New Strategies. Wiley.\cr
Lombardo R Beh EJ 2016 Variants of Simple Correspondence Analysis. The R Journal, 8 (2), 167--184.\cr
Lombardo R Beh EJ and Kroonenberg PM 2016 Modelling Trends in Ordered Correspondence Analysis Using Orthogonal
Polynomials. Psychometrika, 81(2), 325--349.
}
\author{
Rosaria Lombardo and Eric J. Beh
}
\note{
This function is called from the main plot function \code{plot.CAvariants}.
}
\keyword{ graph }% __ONLY ONE__ keyword per line
