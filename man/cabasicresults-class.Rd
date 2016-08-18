\name{cabasicresults-class}
\Rdversion{1.1}
\docType{class}
\alias{cabasicresults-class}
\alias{cabasicresults}
\title{The class of common objects }
\description{
The class of objects common to the six variants of correspondence analysis.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{cabasicresults(...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{RX}}{Object of class \code{"matrix"}. Row profile in symmetrical variants, transposed of column profile in non symmetrical variants of correspondence analysis}
    \item{\code{CX}}{Object of class \code{"matrix"}. Column profile in symmetrical and non symmetrical variants of correspondence analysis}
    \item{\code{Rweights}}{Object of class \code{"matrix"}. Row weights}
    \item{\code{Cweights}}{Object of class \code{"matrix"}. Column weights}
    \item{\code{Raxes}}{Object of class \code{"matrix"}. Row axes}
    \item{\code{Caxes}}{Object of class \code{"matrix"}. Column axes}
    \item{\code{mu}}{Object of class \code{"numeric"}. Singular values of row space}
    \item{\code{mu2}}{Object of class \code{"numeric"}. Singular values of column space}
    \item{\code{catype}}{Object of class \code{"character"}. Type of variants of correspondence analysis}
    \item{\code{tauDen}}{Object of class \code{"numeric"}. Denominator of the tau index in non symmetrical variants of correspondence analysis}
    \item{\code{Z}}{Object of class \code{"matrix"}. Generalised correlation matrix in ordered variants of symmetrical and non symmetrical correspondence analysis}
    \item{\code{ZtZ}}{Object of class \code{"matrix"}. Squared generalised correlation in the row space}
    \item{\code{tZZ}}{Object of class \code{"matrix"}. Squared generalised correlation in the column space}
  }
}
%\section{Extends}{
%Class \code{"\linkS4class{cacorporateplus}"}, directly.
%}
\section{Methods}{
No methods defined with class \code{cabasicresults} in the signature.
}
\references{
Beh EJ and Lombardo R 2014 Correspondence analysis, Theory, Practice and New Strategies. John Wiley & Sons.
}
\author{
Rosaria Lombardo and Eric J. Beh
}
