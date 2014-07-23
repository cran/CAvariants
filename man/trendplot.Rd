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
The parameter for setting the size of character labels of points in graphical displays.
}
  \item{cex.lab}{
The parameter for setting the size of character labels of axes in graphical displays.
}
  \item{main}{
The title of the graphical display.}
  \item{prop}{
The scaling parameter for specifying the limits of the plotting area. By default, it is equal to 1. 
}
  \item{posleg}{
The parameter for specifying the position of the legend in the graphical function \code{trendplot}. \cr
By default, \code{pos = "topleft"}.}
  \item{xlab}{
The parameter for setting the  character label of the horizontal axis in graphical displays.
}
  \item{ylab}{
The parameter for setting the  character label of the vertical axis in graphical displays.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
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
Rosaria Lombardo and Eric J. Beh
}
\note{
This function is called from the main plot function \code{plot.cacorporateplus}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.

%## The function is currently defined as
%function (f, g, cex = 1, cex.lab = 0.8, main = " ", prop = 0.5, 
%    position = "bottomright", xlab = "First Principal Axis", 
%    ylab = "Second Principal Axis") 
%{
%    nrows <- dim(g)[[1]]
%    ncols <- dim(g)[[2]]
%    leg.txt <- dimnames(g)[[1]]
%    colsymb <- c(1:nrows)
%    gt <- t(g)
%    plot(f, g[1, ], type = "b", ylim = range(gt[1:ncols, ], g[1:nrows, 
%        ])/prop, xlab = xlab, ylab = ylab, cex = cex, cex.lab = cex.lab, 
%        main = main, col = 1)
%    for (i in 1:(nrows)) {
%        lines(f, g[i, ], type = "b", pch = i, col = i)
%    }
%    legend(x = position, legend = leg.txt, col = colsymb, pch = c(1:nrows), 
%        bty = "o", cex = 0.8)
%  }
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
\keyword{ graph }% __ONLY ONE__ keyword per line
