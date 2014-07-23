\name{sonscabasic}
\alias{sonscabasic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Singly or one-way ordered non symmetrical correspondence analysis, it implies an ordered (column) variable
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. 
This function is used in the main function \code{CAvariants} when the input parameter is \code{catype="SONSCA"}.
It performs the Hybrid decomposition of the Goodman-Kruskal index and
compute the principal axes for rows and polynomial axes for columns, coordinates, weights of rows and columns, 
inertia as equal to the Pearson's index, the rank of the matrix.
It allows the analyst to decompose the inertia in column polynomial components 
}
\usage{
sonscabasic(xo, mj)
}
\arguments{
\item{xo}{
The two-way contingency table.
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
Beh EJ and Lombardo R 2014 Correspondence analysis, Theory, Practice and New Strategies. John Wiley & Sons.
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
asbestos<-matrix(c(310, 36, 0, 0, 212, 158, 9, 0, 21, 35, 17, 4, 25, 102,  
49, 18, 7, 35, 51, 28), 4, 5, dimnames = list(c("none","grade1", "grade2", "grade3"), 
c("0-9", "10-19", "20-29", "30-39", "40")))
mj<-c(1,2,3,4,5)
sonscabasic(asbestos,mj)
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%
%## The function is currently defined as
%function (xo, mj) 
%{
%    rmax <- min(dim(xo)) - 1
%    x <- xo/sum(xo)
%    rsums <- as.matrix(rowSums(x))
%    csums <- as.vector(colSums(x))
%    tauden <- 1 - sum(rsums^2)
%    drm1 <- diag(1/(rsums + (rsums == 0)) * (1 - (rsums == 0)))
%    dcm1 <- diag(1/(csums + (csums == 0)) * (1 - (csums == 0)))
%    drmh <- diag(rep(1, nrow(x)))
%    dcmh <- sqrt(dcm1)
%    dj <- diag(csums)
%    di <- diag(rsums)
%    uni <- matrix(1, 1, ncol(x))
%    uni1 <- rep(1, nrow(x))
%    Bpoly <- emerson.poly(mj, csums)
%    Bpoly2 <- sqrt(dj) \%*\% Bpoly
%    pcc <- 1/sqrt(tauden) * (x \%*\% dcm1 - rsums \%*\% (uni))
%    u <- svd(pcc)$u
%    Z <- t(u) \%*\% pcc \%*\% dj \%*\% Bpoly
%    ZtZ <- Z \%*\% t(Z)
%    tZZ <- t(Z) \%*\% Z
%    mu2 <- diag(tZZ)
%    mu <- diag(ZtZ)
%    tau <- sum(mu)
%    r <- sum(mu > 1e-15)
%    if (r < rmax) {
%        mu[(r + 1):rmax] <- 0
%    }
%    sonscabasic <- new("cabasicresults", RX = pcc, CX = t(pcc), 
%        Rweights = dj, Cweights = diag(uni1), Raxes = Bpoly, 
%        Caxes = u, r = r, mu = mu, mu2 = mu2, tau = tau, tauDen = tauden, 
%        catype = "SONSCA", Z = Z, ZtZ = ZtZ, tZZ = tZZ)
%  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~multivariate }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
