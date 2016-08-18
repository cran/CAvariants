\name{nscabasic}
\alias{nscabasic}
\title{
Non symmetrical two-way correspondence analysis
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function is used in the main function \code{CAvariants} when the input parameter is \code{catype = "NSCA"}. 
It calculates the singular value decomposition of  the numerator of the Goodman-Kruskal tau index (index of predictability), 
computes principal axes, coordinates, weights of rows and columns, total inertia (equal to the numerator of the tau index) and the rank of the matrix. 
}
\usage{
nscabasic(Xtable)
}
\arguments{
\item{Xtable}{
The two-way contingency table.
}
}
\references{
Beh EJ and Lombardo R 2014 Correspondence Analysis: Theory, Practice and New Strategies. John Wiley & Sons.
}
\author{
Rosaria Lombardo and Eric J. Beh
}
\note{
This function belongs to the \code{R} object class called  \code{cabasicresults}.
}
\examples{
data(asbestos)
nscabasic(asbestos)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate} % __ONLY ONE__ keyword per line
