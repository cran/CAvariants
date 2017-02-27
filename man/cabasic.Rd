\name{cabasic}
\alias{cabasic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Classical two-way correspondence analysis
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function is used in the main function \code{CAvariants} when the input parameter is \code{catype = "CA"}. 
It performs the singular value decomposition of the Pearson's ratio and
compute principal axes, coordinates, weights of rows and columns, total inertia (equal to the Pearson's index) and
the rank of the matrix. 
}
\usage{
cabasic(Xtable)
}
\arguments{
\item{Xtable}{
The two-way contingency table.
}
}
\references{
Beh EJ and Lombardo R 2014 Correspondence Analysis: Theory, Practice and New Strategies. John Wiley & Sons.
Lombardo R Beh EJ Variants of Simple Correspondence Analysis. The R Journal (accepted).
}
\author{
Rosaria Lombardo and Eric J. Beh
}
\note{
This function belongs to the \code{R} object class called  \code{cabasicresults}.
}
\examples{
data(asbestos)
cabasic(asbestos)
}
\keyword{ multivariate}
