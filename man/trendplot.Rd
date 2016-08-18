\name{trendplot}
\alias{trendplot}
\title{
Trends of matrix rows and columns 
}
\description{
This function portrays the row and column trends of the centred column profile matrix reconstructed
by means of othogonal polynomials and/or principal axes.}
\usage{
trendplot(f, g, cex = 1, cex.lab = 0.8, main = " ", prop = 0.5, 
posleg = "topleft", xlab = "First Principal Axis", 
ylab = "Second Principal Axis")
}
\arguments{
  \item{f}{
The row coordinates.
}
  \item{g}{
The column coordinates.
}
  \item{cex}{
The parameter for setting the size of character labels of points in graphical displays. By default, \code{cex = 1}.
}
  \item{cex.lab}{
The parameter for setting the size of character labels of axes in graphical displays. By default, \code{cex.lab = 0.8}.
}
  \item{main}{
The title of the graphical display.}
  \item{prop}{
The scaling parameter for specifying the limits of the plotting area. By default, \code{prop = 0.5}. 
}
  \item{posleg}{
The parameter for specifying the position of the legend in the graphical function \code{trendplot}. 
By default, \code{pos = "topleft"}.}
  \item{xlab}{
The parameter for setting the  character label of the horizontal axis in graphical displays.
}
  \item{ylab}{
The parameter for setting the  character label of the vertical axis in graphical displays.
}
}
\references{
Beh EJ and Lombardo R 2014 Correspondence Analysis: Theory, Practice and New Strategies. John Wiley & Sons.
}
\author{
Rosaria Lombardo and Eric J. Beh
}
\note{
This function is called from the main plot function \code{plot.CAvariants}.
}
\keyword{ graph }% __ONLY ONE__ keyword per line
