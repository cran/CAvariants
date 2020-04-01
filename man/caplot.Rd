\name{caplot}
\alias{caplot}
\title{
Row isometric biplot or Column isometric biplot
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function is used in the main plot function  when the plot type parameter is 
\code{plottype = "biplot"}. It can produce a row or a column biplot. 
}
\usage{
caplot(frows, gcols, firstaxis, lastaxis, nseg, inertiapc, thingseg, col1, col2, 
col3, size1, size2)
}
\arguments{
\item{frows}{
The row principal or standard coordinates.
}
\item{gcols}{
The column principal or standard coordinates.
}
\item{firstaxis}{
The first axis number.
}
\item{lastaxis}{
The second axis number. 
}
\item{nseg}{
The vectors/arrows number where to project  principal (or standard) coordinates. 
}
\item{inertiapc}{
The percentage of the explained inertia by each dimension.
}
\item{thingseg}{
The principal or standard coordinates used to draw vectors (arrows). 
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
\keyword{biplot}
