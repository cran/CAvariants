\name{compsonetable.exe}
\alias{compsonetable.exe}
\title{
Polynomial component of inertia in column space
}
\description{
This function allows the analyst to compute the contribution of the polynomial components to the inertia (chi-squared or tau).
The ordered variable should be the column variable that is transformed by polynomials.
The polynomial components are the column polynomial components. 
The given input matrix is the Z matrix of generalised correlations from the Hybrid decomposition.
It is called by \code{CAvariants} when \code{catype="SOCA"} or \code{"SONSCA"}
}
\usage{
compsonetable.exe(Z)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{Z}{
The matrix of generalised correlations between the polynomial and principal axes.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
The value returned is the matrix
\item{comps}{the matrix of the polynomial component of inertia.}
}
\references{
Beh EJ and Lombardo R 2014 Correspondence analysis, Theory, Methods and New Strategies. Wiley.
}
\author{
Rosaria Lombardo and Eric J. Beh
}
\note{
This function belongs to the class called  \code{cacorporateplus}.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%\examples{
%In case in the main function CAvariants you ask for singly ordered correspondence analysis
%CAvariants((asbestos),catype="SOCA",plottype="classic",scaleplot=1,pos=4)
%The output that this code allows is, for example, with respect to the asbestos data: 
% Polynomial Components of Inertia 
%                        Component Value   P-value
%** Column Components **        0.000000 0.0000000
%Location                     546.399593 0.0000000
%Dispersion                    99.605164 0.0000000
%Error                          2.806743 0.9458938
%** C-Statistic **            648.811501 0.0000000
%## The function is currently defined as
%function (Z) 
%{
%    ZtZ <- Z \%*\% t(Z)
%    factor <- sum(diag(ZtZ))
%    comps <- matrix(0, nrow = 5, ncol = 2)
%    comps[2, 1] <- ZtZ[1, 1]
%    comps[2, 2] <- 1 - pchisq(comps[1, 1], ncol(Z))
%    comps[3, 1] <- ZtZ[2, 2]
%    comps[3, 2] <- 1 - pchisq(comps[2, 1], ncol(Z))
%    comps[4, 1] <- factor - (comps[2, 1] + comps[3, 1])
%    if (nrow(Z) > 2) {
%        comps[4, 2] <- 1 - pchisq(comps[4, 1], (nrow(Z) - 2) * 
%            ncol(Z))
%    }
%    else {
%        comps[4, 2] <- 0
%    }
%    comps[5, 1] <- factor
%    comps[5, 2] <- 1 - pchisq(comps[5, 1], nrow(Z) * ncol(Z))
%    return(comps)
%  }
%}
%% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
