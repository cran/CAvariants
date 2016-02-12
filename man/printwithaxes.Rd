\name{printwithaxes}
\alias{printwithaxes}
\title{
Secondary printing function 
}
\description{
The function is called from the main print function \code{print.cacorporateplus}. 
It adds the names to objects. 
}
\usage{
printwithaxes(x, thenames,digits=3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An \code{R} object.}
  \item{thenames}{
A character vector of up to the same length as \code{x}.}
\item{digits}{ The minimum number of significant digits to be printed in values. Per default, it is equal to 3.}
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
Beh EJ and Lombardo R 2014 Correspondence Analysis: Theory, Practice and New Strategies. John Wiley & Sons.
}
\author{
Rosaria Lombardo and Eric J. Beh
}
\note{
It is called from \code{print.cacorporateplus}.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.

%## The function is currently defined as
%function (res, thenames) 
%{
%    names(res) <- thenames
%    print(res, digits = 4)
%  }
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
