\name{cacorporateplus-class}
\Rdversion{1.1}
\docType{class}
\alias{cacorporateplus-class}
\alias{cacorporateplus}
\title{Class \code{"cacorporateplus"}}
\description{
The Class of objects to show the graphical and numeric results of one of the six variants of correspondence analysis
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{cacorporateplus(...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{S}:}{Object of class \code{"cabasicresults"} The objects of the basic type of correspondence analysis performed}
    \item{\code{DataMatrix}:}{Object of class \code{"matrix"} Data matrix}
    \item{\code{rows}:}{Object of class \code{"numeric"} Number of rows}
    \item{\code{cols}:}{Object of class \code{"numeric"} Number of columns}
    \item{\code{rowlabels}:}{Object of class \code{"character"} Labels of rows}
    \item{\code{collabels}:}{Object of class \code{"character"} Label of columns}
    \item{\code{Rprinccoord}:}{Object of class \code{"matrix"} Row principal coordinates}
    \item{\code{Cprinccoord}:}{Object of class \code{"matrix"} Column principal coordinates}
    \item{\code{Rstdcoord}:}{Object of class \code{"matrix"} Row standard coordinates}
    \item{\code{Cstdcoord}:}{Object of class \code{"matrix"} Column standard coordinates}
    \item{\code{inertiasum}:}{Object of class \code{"numeric"} Total inertia}
    \item{\code{inertias}:}{Object of class \code{"matrix"} Inertia, Inertia percentages, Cumulative inertia percentages of Row space}
    \item{\code{inertias2}:}{Object of class \code{"matrix"} Inertia, Inertia percentages, Cumulative inertia percentages of Column space}
    \item{\code{comps}:}{Object of class \code{"matrix"} Sources of inertia when an ordered analysis is performed}
    \item{\code{printdims}:}{Object of class \code{"numeric"} Number of dimension to print }
    \item{\code{maxaxes}:}{Object of class \code{"numeric"} Maximum number of axes}
    \item{\code{catype}:}{The type of correspondence analysis variant considered}
    \item{\code{mj}:}{Object of class \code{"numeric"} Ordered scores of column variable when different from natural scores}
    \item{\code{mi}:}{Object of class \code{"numeric"} Ordered scores of row variable when different from natural scores}
    \item{\code{pcc}:}{Object of class \code{"matrix"} Centred column profile matrix}
    \item{\code{Jmass}:}{Object of class \code{"matrix"} Column weights}
    \item{\code{Imass}:}{Object of class \code{"matrix"} Row weights}
    \item{\code{Trend}:}{Object of class \code{"matrix"} Trend of rows or columns when an ordered analysis is performed}
    \item{\code{Z}:}{Object of class \code{"matrix"} Matrix of generalised correlations in case of ordinal analysis}
  }
}
\section{Methods}{
No methods defined with class "cacorporateplus" in the signature.
}
\references{
Beh EJ and Lombardo R 2014 Correspondence analysis, Theory, Practice and New Strategies. John Wiley & Sons.
Lombardo R and Beh EJ 2014 Variants of Simple Correspondence Analysis. Submitted to Journal of Statistical Software.

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
%showClass("cacorporateplus")
%}
%\keyword{classes}
%\usage{
%cacorporateplus(...)
%}
%\arguments{
%  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
%}
%}
