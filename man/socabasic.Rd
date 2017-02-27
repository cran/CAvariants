\name{socabasic}
\alias{socabasic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Singly or one-way ordered correspondence analysis, it implies an ordered (column) variable
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function is used in the main function \code{CAvariants} when the input parameter is \code{catype = "SOCA"}.
It performs the hybrid decomposition of the Pearson's ratio and
computes principal axes for rows and polynomial axes for columns, coordinates, weights of rows and columns, 
total inertia (equal to the Pearson's index) and the rank of the matrix.
It allows the analyst to decompose the inertia in terms of the column polynomial components. 
}
\usage{
socabasic(Xtable, mj)
}
\arguments{
\item{Xtable}{
The two-way contingency table.
}
\item{mj}{
The set of ordered column scores.
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
mj <- c(1,2,3,4,5)
socabasic(asbestos, mj)
}
\keyword{ multivariate }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
