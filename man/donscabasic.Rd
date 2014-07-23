\name{donscabasic}
\alias{donscabasic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Doubly or two-way ordered non symmetrical correspondence analysis, it implies two ordered variables
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. 
This function is used in the main function \code{CAvariants} when the input parameter is \code{catype="DONSCA"}. 
It performs the Bivariate moment decomposition of the tau index of predictability. 
It computes polynomial axes, coordinates, weights of rows and columns, total inertia (equal to the Goodman and Kruskal's tau index) and the rank of the matrix. 
It allows the analyst to decompose the inertia into row and column polynomial components. 
}
\usage{
donscabasic(xo, mi, mj)
}
\arguments{
\item{xo}{
The two-way contingency table.
}
\item{mi}{
The set of ordered row scores.
}
\item{mj}{
The set of ordered  column scores.
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
This function belongs to the \code{R} object class called \code{cabasicresults}.
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
donscabasic(asbestos, mi, mj)
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%
%## The function is currently defined as
%function (xo, mi, mj) 
%{
%    rmax <- min(dim(xo)) - 1
%    x <- xo/sum(xo)
%    rsums <- as.vector(rowSums(x))
%    csums <- as.vector(colSums(x))
%    drm1 <- diag(1/(rsums + (rsums == 0)) * (1 - (rsums == 0)))
%    dcm1 <- diag(1/(csums + (csums == 0)) * (1 - (csums == 0)))
%    drmh <- diag(rep(1, nrow(x)))
%    dcmh <- sqrt(dcm1)
%    dj <- diag(csums)
%    di <- diag(rsums)
%    tauden <- 1 - sum(rsums^2)
%    Apoly <- emerson.poly(mi, rsums)
%    Apoly2 <- sqrt(di) \%*\% Apoly
%    Bpoly <- emerson.poly(mj, csums)
%    Bpoly2 <- sqrt(dj) \%*\% Bpoly
%    pcc <- 1/sqrt(tauden) * (drmh \%*\% (x - rsums \%*\% t(csums)) \%*\% 
%        dcm1)
%    Z <- t(Apoly2) \%*\% pcc \%*\% dj \%*\% Bpoly
%    pi <- (Apoly2) \%*\% Z \%*\% t(Bpoly2)
%    ZtZ <- Z \%*\% t(Z)
%    mu <- diag(ZtZ)
%    tau <- sum(mu)
%    tZZ <- t(Z) \%*\% Z
%    mu2 <- diag(tZZ)
%    r <- sum(mu > 1e-15)
%    Cweights <- dj
%    donscabasic <- new("cabasicresults", RX = pcc, CX = t(pcc), 
%        Rweights = dj, Cweights = diag(rep(1, nrow(xo))), Raxes = Bpoly, 
%        Caxes = Apoly2, r = r, mu = mu, mu2 = mu2, catype = "DONSCA", 
%        tau = tau, tauDen = tauden, Z = Z, ZtZ = ZtZ, tZZ = tZZ)
%  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~multivariate }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
