\name{emerson.poly}
\alias{emerson.poly}
\title{
Orthogonal polynomials}
\description{
This function is called  from the functions
\code{docabasic, socabasic, sonscabasic} and \code{donscabasic}. 
It allows the analyst to compute the orthogonal polynomials of the ordered categorical variable. 
The number of the polynomials is equal to the variable category less one.
The function computes the polynomial transformation of the ordered categorical variable.}
\usage{
emerson.poly(mj, pj)
}
\arguments{
  \item{mj}{
 The ordered scores of an ordered variable. By default \code{mj=NULL}, the natural scores (1,2,...) are computed.
}
  \item{pj}{
 The marginals, relative frequencies of the ordered variable. 
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Describe the value returned
\item{B}{the matrix of the orthogonal polynomials without the trivial polynomial.}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Beh EJ and Lombardo R 2014 Correspondence analysis, Theory, Methods and New Strategies. Wiley.
}
\author{
Rosaria Lombardo and Eric J Beh
}
\note{
Note that the sum of the marginals of the ordered variables should be one.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
emerson.poly(c(1,2,3,4,5),as.vector(c(.1,.2,.3,.2,.2)))
%           [,1]       [,2]        [,3]       [,4]
%[1,] -1.7614097  1.9490178 -1.32656790  0.5822225
%[2,] -0.9607689 -0.1461763  1.30369604 -1.1644450
%[3,] -0.1601282 -0.9745089  0.04574372  1.1644450
%[4,]  0.6405126 -0.5359799 -1.39518348 -1.1644450
%[5,]  1.4411534  1.1694107  0.68615581  0.2911113
%## The function is currently defined as
%function (mj, pj) 
%{
%    nc <- length(mj)
%    Dj <- diag(pj)
%    B <- matrix(1, (nc + 1), nc)
%    B[1, ] <- 0
%    Sh <- Th <- Vh <- NULL
%    for (i in 3:(nc + 1)) {
%        for (j in 1:nc) {
%            Th[i] <- mj \%*\% Dj \%*\% B[i - 1, ]^2
%            Vh[i] <- mj \%*\% Dj \%*\% (B[i - 1, ] * B[i - 2, ])
%            Sh[i] <- sqrt(mj^2 \%*\% Dj \%*\% B[i - 1, ]^2 - Th[i]^2 - 
%                Vh[i]^2)^(-1)
%            B[i, j] <- Sh[i] * ((mj[j] - Th[i]) * B[i - 1, j] - 
%                Vh[i] * B[i - 2, j])
%        }
%    }
%    B <- t(B)
%    B <- B[, -c(1, 2)]
%    return(B)
%  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
\keyword{ nonparametric }% __ONLY ONE__ keyword per line
