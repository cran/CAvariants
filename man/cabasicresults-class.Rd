\name{cabasicresults-class}
\Rdversion{1.1}
\docType{class}
\alias{cabasicresults-class}
\alias{cabasicresults}

\title{Class \code{"cabasicresults"}}
\description{
The Class of objects common to the six variants of correspondence analysis
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{cabasicresults(...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{RX}:}{Object of class \code{"matrix"} Row profile in symmetrical variants, transposed of column profile in non symmetrical variants of correspondence analysis}
    \item{\code{CX}:}{Object of class \code{"matrix"} Column profile in symmetrical and non symmetrical variants of correspondence analysis}
    \item{\code{Rweights}:}{Object of class \code{"matrix"} Row weights}
    \item{\code{Cweights}:}{Object of class \code{"matrix"} Column weights}
    \item{\code{Raxes}:}{Object of class \code{"matrix"} Row axes}
    \item{\code{Caxes}:}{Object of class \code{"matrix"} Column axes}
 %   \item{\code{r}:}{Object of class \code{"numeric"} Rank of the matrix}
    \item{\code{mu}:}{Object of class \code{"numeric"} Singular values of row space}
    \item{\code{mu2}:}{Object of class \code{"numeric"} Singular values of column space}
    \item{\code{catype}:}{Object of class \code{"character"} Type of variants of correspondence analysis}
 %   \item{\code{tau}:}{Object of class \code{"numeric"} Tau index in non symmetrical variants of correspondence analysis }
    \item{\code{tauDen}:}{Object of class \code{"numeric"} Denominator of tau index in non symmetrical variants of correspondence analysis}
    \item{\code{Z}:}{Object of class \code{"matrix"} Generalised correlation matrix in ordered variants of symmetrica and non symmetrical correspondence analysis}
    \item{\code{ZtZ}:}{Object of class \code{"matrix"} Squared generalised correlation in Row space}
    \item{\code{tZZ}:}{Object of class \code{"matrix"} Squared generalised correlation in Column space}
    %\item{\code{S}:}{Object of class \code{"cabasicresults"} List of objects of Class cabasicresults}
    %\item{\code{DataMatrix}:}{Object of class \code{"matrix"} Data table}
    %\item{\code{rows}:}{Object of class \code{"numeric"} Number of rows}
    %\item{\code{cols}:}{Object of class \code{"numeric"} Number of columns}
    %\item{\code{rowlabels}:}{Object of class \code{"character"} Labels of rows}
    %\item{\code{collabels}:}{Object of class \code{"character"} Labels of columns}
    %\item{\code{Rprinccoord}:}{Object of class \code{"matrix"} Row principal coordinates}
    %\item{\code{Cprinccoord}:}{Object of class \code{"matrix"} Column principal coordinates}
    %\item{\code{Rstdcoord}:}{Object of class \code{"matrix"} Row standard coordinates}
    %\item{\code{Cstdcoord}:}{Object of class \code{"matrix"} Column standard coordinates}
    %\item{\code{inertiasum}:}{Object of class \code{"numeric"} Total inertia}
    %\item{\code{inertias}:}{Object of class \code{"matrix"} Inertia, inertia percentage, cumulative inertia percentage in Row space}
    %\item{\code{inertias2}:}{Object of class \code{"matrix"} Inertia, inertia percentage, cumulative inertia percentage in Column space }
    %\item{\code{comps}:}{Object of class \code{"matrix"} Sources of inertia when an ordered analysis is performed}
    %\item{\code{printdims}:}{Object of class \code{"numeric"} Number of dimension to print}
    %\item{\code{maxaxes}:}{Object of class \code{"numeric"} Maximum number of axes}
    %\item{\code{mj}:}{Object of class \code{"numeric"} Ordered scores of column variable when different from natural scores}
    %\item{\code{mi}:}{Object of class \code{"numeric"} Ordered scores of row variable when different from natural scores}
    %\item{\code{pcc}:}{Object of class \code{"matrix"} Centred column profile matrix}
    %\item{\code{Jmass}:}{Object of class \code{"matrix"} Column weights}
    %\item{\code{Imass}:}{Object of class \code{"matrix"} Row weights}
    %\item{\code{Trend}:}{Object of class \code{"matrix"} Trend of rows or columns when an ordered analysis is performed}
  }
}
%\section{Extends}{
%Class \code{"\linkS4class{cacorporateplus}"}, directly.
%}
\section{Methods}{
No methods defined with class "cabasicresults" in the signature.
}
\references{
Beh EJ and Lombardo R 2014 Correspondence analysis, Theory, Practice and New Strategies. John Wiley & Sons.
}
\author{
Rosaria Lombardo and Eric J. Beh
}
%\note{
%%  ~~further notes~~
%}
%% ~Make other sections like Warning with \section{Warning }{....} ~
%\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
%}
%\examples{
%showClass("cabasicresults")
%}
%\keyword{classes}
%\usage{
%cabasicresults(...)
%}
%\arguments{
%  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
%}
%}
