\name{caplotord}
\alias{caplotord}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Row isometric or column isometric biplot for ordered variants of correspondence analysis
}
\description{
This function is used in the main plot function when the plot type parameter is 
\code{plottype = "biplot"}. It can produce a row polynomial biplot or a column polynomial biplot. 
}
\usage{
caplotord(frows, gcols, firstaxis, lastaxis, nseg, inertiapc, thingseg, col1, 
col2, col3, size1, size2)
}
\arguments{
\item{frows}{
The row principal or standard coordinates.
}
\item{gcols}{
The column principal or standard coordinates.
}
\item{firstaxis}{
The first polynomial axis number. 
}
\item{lastaxis}{
The second polynomial axis number. 
}
\item{nseg}{
The vectors/arrows number where to project  principal (or standard) coordinates. 
}
\item{inertiapc}{
The percentage of the explained inertia by each dimension.
}
\item{thingseg}{
The   principal or standard coordinates used to draw vectors (arrows). 
}
\item{col1}{
The colour for the row variable labels.
}
\item{col2}{
The colour for the column variable labels.
}
\item{col3}{
The colour for the  vectors (arrows) used in biplots.
}
\item{size1}{
The size of the plotted symbol for categories in biplot. 
}
\item{size2}{
The size of the plotted text for categories in biplot. 
}
}
\references{
Beh EJ and Lombardo R 2014 Correspondence Analysis: Theory, Practice and New Strategies. Wiley.\cr
Lombardo R Beh EJ 2016 Variants of Simple Correspondence Analysis. The R Journal, 8 (2), 167--184.
}
\author{
Rosaria Lombardo and Eric J. Beh
}
\note{
This function depends on the \code{R} library \code{plotly}.
}
%\examples{
%data(asbestos)
%cabasic(asbestos)
%}
\keyword{multivariate}
