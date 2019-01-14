\name{caplot3d}
\alias{caplot3d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Three dimensional correspondence plot
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function is used in the plot function \code{plot.CAvariants} when the logical parameter is 
\code{plot3d = TRUE}. 
It produces a 3-dimensional visualization of the association. 
}
\usage{
caplot3d(f, g, percIn, size=10, adj=c(0,0.5) ,prop=1)
}
\arguments{
\item{f}{
The row principal or standard coordinates.
}
\item{g}{
The column principal or standard coordinates.
}
\item{percIn}{
The percentage of the total inertia explained inertia by each principal axis.
}
\item{size}{
The size of the plotting area. By default, \code{size = 10}.
}
\item{adj}{
One or two values in the interval [0, 1] that specifies the x (and optionally y) position of the label categories. 
By default, \code{adj = c(0, 0.5).}
}
\item{prop}{
The scaling parameter for specifying the limits of the plotting area. By default, \code{prop = 1}.}
}
\references{
Beh EJ and Lombardo R 2014 Correspondence Analysis: Theory, Practice and New Strategies. Wiley.\cr
Lombardo R Beh EJ 2016 Variants of Simple Correspondence Analysis. The R Journal, 8 (2), 167--184.
}
\author{
Rosaria Lombardo and Eric J. Beh
}
\note{
This function depends on the \code{R} library \code{rgl}.
}
%\examples{
%data(asbestos)
%cabasic(asbestos)
%}
\keyword{ multivariate}
