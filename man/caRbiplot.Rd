\name{caRbiplot}
\alias{caRbiplot}
\title{
Row isometric biplot or Column isometric biplot
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function is used in the main plot function  when the plot type parameter is 
\code{plottype = "biplot"}. It can produce a row biplot or a column biplot. 
}
\usage{
caRbiplot(frows, gcols, firstaxis, lastaxis, inertiapc,  bip="row", size1,size2)
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
\item{inertiapc}{
The percentage of the explained inertia by each dimension.
}
\item{bip}{
The type of  biplot. One may specify a row-isometric biplot or a column-isometric biplot (when using 
in  the function \code{plot.CAvariant} the parameter \code{biptype = "row"} or \code{biptype = "column"}). 
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
%data(asbestos)
%cabasic(asbestos)
%}
\keyword{biplot}
