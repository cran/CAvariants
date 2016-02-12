\name{nscabasic}
\alias{nscabasic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Non-symmetrical two-way correspondence analysis
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function is used in the main function \code{CAvariants} when the input parameter is \code{catype="NSCA"}. 
It calculates the singular value decomposition of the index of predictability (Goodman and Kruskal tau index), 
computes principal axes, coordinates, weights of rows and columns, total inertia (equal to the tau index) and the rank of the matrix. 
}
\usage{
nscabasic(X)
}
\arguments{
\item{X}{
The two-way contingency table.
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
nscabasic(asbestos)
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%## The function is currently defined as
%function (X) 
%{
%    X <- X/sum(X)
%    r <- min(dim(X)) - 1
%    Imass <- rowSums(X)
%    tauden <- 1 - sum(Imass^2)
%    rsums <- as.vector(rowSums(X))
%    csums <- as.vector(colSums(X))
%    drm1 <- diag(1/(rsums + (rsums == 0)) * (1 - (rsums == 0)))
%    dcm1 <- diag(1/(csums + (csums == 0)) * (1 - (csums == 0)))
%    drmh <- diag(rep(1, nrow(X)))
%    dcmh <- sqrt(dcm1)
%    Z <- 1/sqrt(tauden) * (drmh \%*\% (X - rsums \%*\% t(csums)) \%*\% 
%        dcmh)
%    Y <- svd(Z)
%    mu <- Y$d
%    tau <- sum(mu^2)
%    R <- drm1 \%*\% X
%    C <- dcm1 \%*\% t(X)
%    nscabasic <- new("cabasicresults", RX = R, CX = C, Rweights = dcmh, 
%        Cweights = drmh, Raxes = dcmh \%*\% Y$v, Caxes = Y$u, r = r, 
%        mu = mu, tau = tau, tauDen = tauden, catype = "NSCA")
%  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
