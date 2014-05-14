\name{compstable.exe}
\alias{compstable.exe}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Polynomial component of inertia in row and column spaces
}
\description{
This function allows the analyst to compute the contribution of the polynomial components to the inertia (chi-squared or tau).
The ordered variable should be the both row and column variables that are transformed by polynomials.
The polynomial components are the row and column polynomial components. 
The given input matrix is the Z matrix of generalised correlations from the Bivariate moment decomposition.
It is called by \code{CAvariants} when \code{catype="DOCA"} or \code{"DONSCA"}
}
\usage{
compstable.exe(Z)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{The matrix of generalised correlations between the polynomial axes.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
The value returned is the matrix
\item{comps }{The matrix of the polynomial component of inertia.}
}
\references{
Beh EJ and Lombardo R 2014 Correspondence analysis, Theory, Methods and New Strategies. Wiley.
}
\author{
Rosaria Lombardo and Eric J. Beh
}
\note{
This function belongs to the \code{R} object class called  \code{cacorporateplus}.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%\examples{
%In case in the main function CAvariants you ask for singly ordered correspondence analysis
%CAvariants((asbestos),catype="DOCA",plottype="classic",scaleplot=1,pos=4)
%The output that this code allows is, for example, with respect to the asbestos data: 
% Polynomial Components of Inertia 
%                            Component Value     P-value
%** Row Components **                0.00000 0.000000000
%Location                          785.59362 0.000000000
%Dispersion                        110.71396 0.000000000
%Error                              21.77819 0.001328172
%** Column Components **             0.00000 0.000000000
%Location                          759.80892 0.000000000
%Dispersion                        142.36277 0.000000000
%Error                              15.91408 0.003136643
%** Chi-squared Statistic **       918.08577 0.000000000
%## The function is currently defined as
%function (Z) 
%{
%    tZZ <- t(Z) \%*\% Z
%    ZtZ <- Z \%*\% t(Z)
%    factor <- sum(diag(tZZ))
%    comps <- matrix(0, nrow = 9, ncol = 2)
%    comps[2, 1] <- tZZ[1, 1]
%    comps[2, 2] <- 1 - pchisq(comps[2, 1], nrow(Z))
%    comps[3, 1] <- tZZ[2, 2]
%    comps[3, 2] <- 1 - pchisq(comps[3, 1], nrow(Z))
%    comps[4, 1] <- factor - (comps[2, 1] + comps[3, 1])
%    if (ncol(Z) > 2) {
%        comps[4, 2] <- 1 - pchisq(comps[4, 1], nrow(Z) * (ncol(Z) - 
%            2))
%    }
%    else {
%        comps[4, 2] <- 0
%    }
%    comps[6, 1] <- ZtZ[1, 1]
%    comps[6, 2] <- 1 - pchisq(comps[6, 1], ncol(Z))
%    comps[7, 1] <- ZtZ[2, 2]
%    comps[7, 2] <- 1 - pchisq(comps[7, 1], ncol(Z))
%    comps[8, 1] <- factor - (comps[6, 1] + comps[7, 1])
%    if (nrow(Z) > 2) {
%        comps[8, 2] <- 1 - pchisq(comps[8, 1], (nrow(Z) - 2) * 
%            ncol(Z))
%    }
%    else {
%        comps[8, 2] <- 0
%    }
%    comps[9, 1] <- factor
%    comps[9, 2] <- 1 - pchisq(comps[9, 1], nrow(Z) * ncol(Z))
%    return(comps)
%  }
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
