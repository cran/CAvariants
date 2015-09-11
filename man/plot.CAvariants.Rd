\name{plot.CAvariants}
\alias{plot.CAvariants}
\title{
Main plot function 
}
\description{
This function allows the analyst to produce the suitable graphical displays with respect to the six variants of correspondence analysis. 
The main plot function called from the main function \code{CAvariants}. 
In particular when \code{plottype = "classic"}, it produces classical graphical displays for \code{catype = "CA"} and \code{catype = "NSCA"}, 
where the row and column variables are graphically depicted in principal coordinates.\cr
When we set \code{plottype = "biplot"}, it produces biplot graphical displays, or polynomial biplots in case of ordered analysis.
Notice that for ordered analysis only polynomial biplots are suitable.  
When the input parameter \code{catype}  is equal to \code{catype = "SOCA"}, 
\code{catype = "DOCA"}, \code{catype = "SONSCA"} or \code{catype = "DONSCA"}, then the input parameter  
\code{plottype} should be equal to \code{plottype="biplot"},  if \code{biptype="row"} , it will give back a row isometric polynomial biplot.
 }
\usage{
\S3method{plot}{CAvariants}(x, firstaxis = 1, lastaxis = 2, cex = 1, 
prop = 1,  plottype = "biplot", biptype = "rows", M = 2, 
scaleplot = 3, posleg = "topleft", pos = 2, ell = FALSE, \dots)
}
\arguments{
  \item{x}{
Represents the set of the output parameters of the main function \code{CAvariants} of the \code{R} object class \code{cacorporateplus}.
}
  \item{firstaxis}{
The dimension reflected along  the horizontal axis. By default, it is equal to 1.
}
  \item{lastaxis}{
The dimension reflected along the vertical axis. By default, it is equal to 2.
}
  \item{cex}{
The parameter that specifies the size of character labels of points in graphical displays. By default, it is equal to 1.
}
  \item{prop}{
The scaling parameter for specifying the limits of the plotting area. By default, it is equal to 1. 
}
  
  \item{plottype}{
 The input parameter specifying what kind of graphical display is requested. By default, it is equal 
to \code{plottype = "biplot"}.
}
\item{biptype}{
The input parameter specifying what kind of biplot is requested. By default, it is equal to \code{biptype = "row"}
(row isometric biplot).
}
\item{M}{
The number of  axes used when drawing ellypse. By default, it is equal to \code{M= 2}.
Its maximum value is equal to the rank of the data matrix.
}
\item{scaleplot}{
The scaling parameter for scaling biplot coordinates. By default, \code{scaleplot = 1}.
}
\item{posleg}{
The parameter for specifying the position of the legend in the graphical function \code{trendplot}. 
By default, \code{posleg= "topleft"}.}
 \item{pos}{
The parameter for specifying the position of  point symbols in the graphical displays \cr
By default, \code{pos=2}.}
\item{ell}{ The flag parameter to allow the building of algebraic ellipses. By default, \code{ell = F}.}
\item{\dots}{Further arguments passed to or from other methods.}
}
\details{
It produces classical and biplot graphical displays. Further when \code{catype} is equal to \cr
\code{"DOCA", "SOCA", "DONSCA"} or \code{"SONSCA"}, 
the trend of row and column variables after the reconstruction of column profiles by polynomials is portrayed.\cr
For classical biplot displays, it superimposes on it algebraic ellipses of confidence. It uses the secondary plot functions \code{caellipse} or
\code{nsca.ellipse}, depending on the input parameter \code{catype}. }
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
\references{
Beh EJ and Lombardo R 2014 Correspondence analysis, Theory, Practice and New Strategies. John Wiley & Sons.
}
\author{
Rosaria Lombardo and Eric J Beh
}
\note{
For classical graphical displays, both sets of coordinates are defined using principal coordinates (see Greenacre's terminology).
In biplot graphical displays, one set of coordinates is standard and the other is principal. 
When the analysis is ordered, it makes sense only biplot. One set of coordinates consists of standard polynomial coordinates and the other one is of principal polynomial coordinates.

}
%% ~Make other sections like Warning with \section{Warning }{....} ~
%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
asbestos<-matrix(c(310, 36, 0, 0, 212, 158, 9, 0, 21, 35, 17, 4, 25,102,  
49, 18, 7, 35, 51, 28), 4, 5, dimnames = list(c("none","grade1", "grade2", "grade3"), 
c("0-9", "10-19", "20-29", "30-39", "40")))
risasbestos<-CAvariants(asbestos,catype = "DOCA",firstaxis=1,lastaxis=2) 
plot.CAvariants(risasbestos,plotype="biplot",biptype="row")
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%## The function is currently defined as
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
\keyword{multivariate}% __ONLY ONE__ keyword per line
