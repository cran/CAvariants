\name{plotone}
\alias{plotone}
\title{
Secondary plot function
}
\description{
This function produces a classical plot when the input parameter of the main code \code{CAvariants} 
has been set as \code{plottype = "classic"}, or it gives a biplot when  \code{plottype = "biplot"}.
This function allows the analyst to plot row and column points in graphical displays, using symbols and labels for each point. 
It is called internally from the main graphical function \code{plot.CAvariants}.
}
\usage{
plotone(a1, a2, plottype, things, nthings, nvars, Thingcoord, Varcoord, 
inertiapc, thinggroup, thinggrlab, vargroup, vargrlab, thinglabels, varlabels, 
picsize, cex = 0.8, cex.lab = 0.8, type , catype , pos = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a1}{
The dimension reflected along the horizontal (polynomial or principal) axis in the graphical display.
}
  \item{a2}{
 The dimension reflected along the vertical (polynomial or principal) axis in the graphical display.
}
  \item{plottype}{
The paramater for specifying what kind of graphical display is required. 
By default, \code{plottype = "classic"}, otherwise \code{plottype = "biplot"}. 
}
  \item{things}{
The input parameter for specifying what variable is in principal coordinates. Its value is \code{"row"} or \code{"column"} in biplots 
\code{"rows&columns"} in classical plots.}
  \item{nthings}{
The number of row categories.}
  \item{nvars}{
The number of  column categories.}
  \item{Thingcoord}{
The coordinates of the row variable.}
  \item{Varcoord}{
The coordinates of the column variable.}
  \item{inertiapc}{
The percentage of explained inertia for each axis.
}
  \item{thinggroup}{
It is a list of two objects, the first is the vector of the ordered natural scores of the row variable, the second object is a vector of one whose length is equal to the row number.
}
  \item{thinggrlab}{
The set of plotting symbols and colours of the row variable.
}
  \item{vargroup}{
It is a list of two objects, the first is the vector of the ordered natural scores of the column variable, the second object is a vector of one whose length is equal to the column number.
}
  \item{vargrlab}{
The set of plotting symbols and colours of the column variable.
}
  \item{thinglabels}{
The labels of the row variable.}
  \item{varlabels}{
The labels of the column variable.}
  \item{picsize}{
The graphical parameter for specifying the size of the plotting area.}
  \item{cex}{
The parameter for setting the size of character labels of points in graphical displays. By default, \code{cex = 0.8}.
}
  \item{cex.lab}{
The parameter that specifies the size of character labels of axes in graphical displays. By default, \code{cex.lab = 0.8}.
}
  \item{type}{
The graphical parameter that specifies the type of line to be drawn.}
  \item{catype}{
The input parameter specifying what variant of correspondence analysis is to be performed.
}
  \item{pos}{
The parameter for specifying the position of  point symbols in the graphical displays. 
By default, \code{pos = 2}.}
}
\references{
Beh EJ and Lombardo R 2014 Correspondence Analysis: Theory, Practice and New Strategies. John Wiley & Sons.
Lombardo R Beh EJ Variants of Simple Correspondence Analysis. The R Journal (accepted).
}
\author{
Rosaria Lombardo and Eric J Beh
}
\note{
In classical graphical displays, both sets of coordinates are principal coordinates (see Greenacre's terminology).
In biplot graphical displays, one set of coordinates is standard and the other is principal. 
When the analysis is ordered, one set of coordinates consists of standard polynomial coordinates and the other one is of principal polynomial coordinates.
This function is internally called by the main graphical function \code{plot.CAvariants}.
}
\keyword{ graph }% __ONLY ONE__ keyword per line
