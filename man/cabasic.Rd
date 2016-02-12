\name{cabasic}
\alias{cabasic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Classical two-way correspondence analysis
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function is used in the main function \code{CAvariants} when the input parameter is \code{catype="CA"}. 
It performs the singular value decomposition of Pearson's ratio and
compute principal axes, coordinates, weights of rows and columns, total inertia (equal to the Pearson's index), 
the rank of the matrix. 
}
\usage{
cabasic(X)
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
This function belongs to the \code{R} object class called  cabasicresults.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
asbestos<-matrix(c(310, 36, 0, 0, 212, 158, 9, 0, 21, 35, 17, 4, 25,102,  
49, 18, 7, 35, 51, 28), 4, 5, dimnames = list(c("none","grade1", "grade2", "grade3"), 
c("0-9", "10-19", "20-29", "30-39", "40")))
cabasic(asbestos)
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%## The function is currently defined as
%function (X) 
%{
%    n <- sum(X)
%    X <- X/n
%    rmax <- min(dim(X)) - 1
%    rsums <- as.vector(rowSums(X))
%    csums <- as.vector(colSums(X))
%    di <- diag(rsums)
%    dj <- diag(csums)
%    drm1 <- diag(1/(rsums + (rsums == 0)) * (1 - (rsums == 0)))
%    dcm1 <- diag(1/(csums + (csums == 0)) * (1 - (csums == 0)))
%    drmh <- sqrt(drm1)
%    dcmh <- sqrt(dcm1)
%    ratio <- (drmh \%*\% (X - rsums \%*\% t(csums)) \%*\% dcmh) * sqrt(n)
%    ratio2 <- drm1 \%*\% X \%*\% dcm1
%    Y <- svd(ratio)
%    mu <- Y$d
%    r <- sum(mu > 1e-15)
%    R <- drm1 \%*\% X
%    C <- dcm1 \%*\% t(X)
%    if (r < rmax) {
%        mu[(r + 1):rmax] <- 0
%        Raxes[, (r + 1):rmax] <- 0
%        Caxes[, (r + 1):rmax] <- 0
%    }
%    cabasic <- new("cabasicresults", RX = R, CX = C, Rweights = dcmh, 
%        Cweights = drmh, Raxes = Y$v, Caxes = Y$u, r = r, mu = mu, 
%        Z = ratio2, catype = "CA")
%  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
