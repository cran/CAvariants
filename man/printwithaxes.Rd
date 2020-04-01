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
A matrix.}
\item{thenames}{
A character vector of the same length as \code{x}.}
\item{digits}{The number of decimal places used for displaying the numerical summaries of the analysis. 
By default, \code{digits = 3}.}
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
Beh EJ and Lombardo R 2014 Correspondence Analysis: Theory, Practice and New Strategies. Wiley.\cr
Lombardo R Beh EJ 2016 Variants of Simple Correspondence Analysis. The R Journal, 8 (2), 167--184.
}
\author{
Rosaria Lombardo and Eric J. Beh
}
\keyword{print}% __ONLY ONE__ keyword per line
