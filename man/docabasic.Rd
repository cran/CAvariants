\name{docabasic}
\alias{docabasic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Doubly or two-way ordered correspondence analysis, it implies two ordered variables
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function is used in the main function \code{CAvariants} when the input parameter is \code{catype="DOCA"}. 
It performs the Bivariate moment decomposition of the  Pearson's ratio, 
computes polynomial axes, coordinates, weights of rows and columns, total inertia (equal to the Pearson's index), the rank of the matrix. \cr
It allows one to decompose the inertia into row and column polynomial components. 
}
\usage{
docabasic(xo, mi, mj)
}
\arguments{
\item{xo}{
The two-way contingency table.
}
\item{mi}{
The set of ordered row scores.
}
\item{mj}{
The set of ordered column scores.
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
Beh EJ and Lombardo R 2014 Correspondence Analysis: Theory, Practice and New Strategies. John Wiley & Sons.
%Beh EJ 1997 Simple Correspondence analysis of ordinal cross-classifications using orthogonal polynomials. 
%Biometrical Journal, 39, 589 -- 613.
}
\author{
Rosaria Lombardo and Eric J. Beh
}
\note{
This function belongs to the \code{R} object class called  \code{cabasicresults}.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
asbestos<-matrix(c(310, 36, 0, 0, 212, 158, 9, 0, 21, 35, 17, 4, 25,102,  
49, 18, 7, 35, 51, 28), 4, 5, dimnames = list(c("none","grade1", "grade2", "grade3"), 
c("0-9", "10-19", "20-29", "30-39", "40")))
mi <- c(1,2,3,4)
mj <- c(1,2,3,4,5)
docabasic(asbestos, mi, mj)
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%
%## The function is currently defined as
%function (xo, mi, mj) 
%{
%    n <- sum(xo)
%    x <- xo/n
%    rsums <- as.vector(rowSums(x))
%    csums <- as.vector(colSums(x))
%    dj <- diag(csums)
%    di <- diag(rsums)
%    drm1 <- diag(1/(rsums + (rsums == 0)) * (1 - (rsums == 0)))
%    dcm1 <- diag(1/(csums + (csums == 0)) * (1 - (csums == 0)))
%    drmh <- sqrt(drm1)
%    dcmh <- sqrt(dcm1)
%    ratio <- drmh \%*\% (x - rsums \%*\% t(csums)) \%*\% dcmh * sqrt(n)
%    R <- drm1 \%*\% x
%    C <- dcm1 \%*\% t(x)
%    rmax <- min(dim(xo)) - 1
%    Bpoly <- emerson.poly(mj, csums)
%    Bpoly2 <- sqrt(dj) \%*\% Bpoly
%    Apoly <- emerson.poly(mi, rsums)
%    Apoly2 <- (di) \%*\% Apoly
%    Z <- t(Apoly2) \%*\% solve(di) \%*\% (ratio) \%*\% (dj) \%*\% (Bpoly)
%    ZZ <- Z^2
%    pi <- (Apoly) \%*\% Z \%*\% t(Bpoly)
%    ZtZ <- Z \%*\% t(Z)
%    tZZ <- t(Z) \%*\% Z
%    mu <- diag(ZtZ)
%    m2 <- diag(tZZ)
%    r <- sum(mu > 1e-15)
%    docabasic <- new("cabasicresults", RX = R, CX = C, Cweights = drmh, 
%        Rweights = dcmh, Caxes = Apoly, Raxes = Bpoly, r = r, 
%        mu = mu, mu2 = diag(tZZ), catype = "DOCA", Z = Z, ZtZ = ZtZ, 
%        tZZ = tZZ)
%  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~multivariate }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
