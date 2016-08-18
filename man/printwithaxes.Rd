\name{printwithaxes}
\alias{printwithaxes}
\title{
Secondary printing function 
}
\description{
The function is called from the main print function \code{print.CAvariants}. 
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
\item{digits}{ The minimum number of significant digits to be printed in values. By default, \code{digits = 3}.}
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
It is called from \code{print.CAvariants}.}
\keyword{ multivariate }% __ONLY ONE__ keyword per line
