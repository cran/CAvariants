\name{CAvariants}
\alias{CAvariants}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Six variants of correspondence analysis 
}
\description{
It performs \cr
1) simple correspondence analysis \cr
2) doubly ordered correspondence analysis \cr
3) singly ordered correspondence analysis  \cr
4) non symmetrical correspondence analysis \cr
5) doubly ordered non symmetrical correspondence analysis \cr
6) singly ordered non symmetrical correspondence analysis 
}
\usage{
CAvariants(Xtable, mj = NULL, mi = NULL, prop = 1, cex = 0.8, 
printdims = 3, firstaxis = 1, lastaxis = 2, catype = "CA",  
plottype = "biplot", biptype="row", scaleplot = 3, posleg="topleft", 
pos = 2, M = min(nrow(Xtable), ncol(Xtable)) - 1,ell = TRUE)
}
\arguments{
  \item{Xtable}{
The two-way contingency table.
}
  \item{mj}{
The ordered scores of a column variable. By default \code{mj=NULL}, the natural ordered numbers are computed.
}
  \item{mi}{
 The ordered scores of a row variable. By default \code{mi=NULL}, the natural ordered numbers are computed.
}
  \item{prop}{
 The scaling parameter for specifying the limits of the plotting area. By default it is equal to 1. 
}
  \item{cex}{
 The scaling parameter for specifying the size of characters in plot.
}
  \item{printdims}{
 The parameter for specifying the dimension of a matrix in print. By default it is equal to 3.
}
  \item{firstaxis}{
 The number of the corresponding horizontal (polynomial or principal) axis in graphical display.
}
  \item{lastaxis}{
 The number of the corresponding vertical (polynomial or principal) axis in graphical display.
}
  \item{catype}{
The input parameter for specifying what variant of correspondence analysis is considered. \cr
By default  \code{catype = "CA"}, other possible values are: \cr
\code{catype = "SOCA", catype = "DOCA", catype = "NSCA", catype = "SONSCA", catype = "DONSCA"}.
}
  \item{plottype}{
 The paramater for specifying what kind of graphical display is required. By default 
\code{plottype="classic"}, otherwise \code{plottype="biplot"}. 
}
  \item{biptype}{
 The parameter for specifying what kind of biplot is equired. By default \code{biptype="row"}, otherwise \code{biptype="column"}.
}
  \item{scaleplot}{
 The scaling parameter for scaling biplot coordinates. By default, \code{scaleplot = 1}.
}
  \item{posleg}{
The parameter for specifying the position of the legend in the graphical function \code{trendplot}. 
By default, \code{pos= "topleft"}.} The location may also be specified by setting pos to a single
keyword from the list \code{"bottomright", "bottom", 
"bottomleft", "left","top", "topright", "right"} and \code{"center"}.
  \item{pos}{
The parameter for specifying the position of  point symbols in the graphical displays. 
By default, \code{pos=2}.}
  \item{M}{ The rank of the two-way contingency table.}
\item{ell}{ The flag parameter to allow the building of algebraic ellipses. By default, \code{ell=TRUE}.}
}
\details{
  ~~ This function belongs to the object class called "cacorporateplus"  ~~
}
\value{
Description of the output returned
\item{DataMatrix }{The two-way contingency table.}
\item{rows }{The row number of the two-way contingency table.}
\item{columns }{The column number of the two-way contingency table.}
\item{rowlabels }{The label of the row variable.}
\item{columnlabels }{The label of the column variable.}
\item{Rprinccoord }{The row principal  coordinates. 
When the input parameter "catype" is \dQuote{DOCA}, \dQuote{SOCA}, \dQuote{SONSCA} or \dQuote{DONSCA} then they are row principal polynomial coordinates.}
\item{Cprinccoord }{The column principal  coordinates. 
When the input parameter "catype" is \dQuote{DOCA}, \dQuote{SOCA}, \dQuote{SONSCA} or \dQuote{DONSCA} then they are column principal polynomial coordinates.}
\item{Rstdcoord }{The row standard  coordinates. 
When the input parameter "catype" is \dQuote{DOCA}, \dQuote{SOCA}, \dQuote{SONSCA} or \dQuote{DONSCA} then they are row standard polynomial coordinates.}
\item{Cstdcoord }{The column standard  coordinates. 
When the input parameter "catype" is \dQuote{DOCA}, \dQuote{SOCA}, \dQuote{SONSCA} or \dQuote{DONSCA} then they are column standard polynomial coordinates.}
\item{inertiasum}{The total inertia of the classical correspondence analysis, when catype is  \dQuote{CA}, \dQuote{DOCA} or \dQuote{SOCA} (the Pearson's index) 
or the inertia of non symmetrical correspondence analysis,  when catype is  \dQuote{NSCA}, \dQuote{DONSCA} or \dQuote{SONSCA}  
(numerator of the Goodman-Kruskal tau index).}
\item{inertias}{The associated inertia in absolute value and percentage, in the row space for each principal or polynomial axis.}
\item{inertias2}{The associated inertia in absolute value and percentage, in the column space for each principal or polynomial axis. 
When catype is CA or NSCA the associated inertia in the row and column spaces are the same for each principal axis.}
\item{comps}{The polynomial components of inertia when the variables are ordered. The inertia of row and/or column space is partitioned in terms
of polynomial components.}
\item{printdims}{
 The  dimension of a matrix in print. By default it is equal to 3.}
\item{maxaxes}{The maximum number of allowed axes.}
\item{catype}{The kind of correspondence analysis chosen.}
\item{mj}{The ordered scores of a column variable. When mj=NULL, the natural ordered numbers are shown.}
\item{mi}{The ordered scores of a row variable. When mi=NULL, the natural ordered numbers are shown.}
\item{pcc}{The weighted centered column profile matrix.}
\item{Jmass}{The weight matrix of the column variable.}
\item{Imass}{The weight matrix of the row variable.}
\item{Trend}{The inner product matrix of biplot coordinates.}
}
\references{
Beh EJ and Lombardo R 2014 Correspondence analysis, Theory, Practice and New Strategies. John Wiley & Sons. \cr
Lombardo R and Beh EJ 2014 Variants of Simple Correspondence Analysis. Submitted to Journal of Statistical Software.
}
\author{
Rosaria Lombardo and Eric J Beh
}
\note{
This function recall internally many other functions, depending on the setting of the input parameter catype, it recalls
one of the six functions which does a variant of correspndence analysis.
After performing a variant of correpsondence analysis, it recall the code for printing and plotting the results. These two
important functions are print.cacorporateplus and plot.cacorporateplus. 
This function belongs to the class cacorporateplus. 
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
asbestos<-matrix(c(310, 36, 0, 0, 212, 158, 9, 0, 21, 35, 17, 4, 25,102,  
49, 18, 7, 35, 51, 28), 4, 5, dimnames = list(c("none","grade1", "grade2", "grade3"), 
c("0-9", "10-19", "20-29", "30-39", "40")))
CAvariants(asbestos,catype = "CA", plottype = "classic", prop = 0.8, 
scaleplot = 3, ell = TRUE, posleg = "topleft") 
shopdataM<-structure(c(81, 66, 150, 667, 67, 24, 47, 430, 743, 132, 32, 
197, 209, 138, 204, 340, 1409, 259, 272, 117, 637, 684, 408, 
57, 547, 550, 304, 193, 229, 527, 258, 368, 98, 246, 116, 298, 
61, 402, 454, 384, 149, 151, 84, 146, 141, 61, 40, 13, 71, 52, 
138, 252, 942, 297, 313, 92, 251, 167, 193, 30, 16, 130, 111, 
280, 624, 359, 109, 136, 36, 96, 67, 75, 11, 16, 31, 54, 200, 
195, 178, 53, 121, 36, 48, 29, 50, 5, 6, 14, 41, 152, 88, 137, 
68, 171, 37, 56, 27, 55, 17, 3, 11, 50, 211, 90, 45, 28, 145, 
17, 41, 7, 29, 28, 8, 10, 28, 111, 34), .Dim = c(13L, 9L), .Dimnames = list(
c("clothing", "accessories", "tobacco", "stationary", "books", 
"records", "household", "candy", "toys", "jewelry", "perfumes", 
"hobby", "other"), c("M12<", "M13", "M16", "M19", "M25", 
"M35", "M45", "M57", "M65+")))
CAvariants(shopdataM, catype = "CA", plottype = "biplot",biptype="row", scaleplot = 4,ell=TRUE)
%%## The function is currently defined as
%function (Xtable, mj = NULL, mi = NULL, prop = 1, cex = 0.8, 
%    printdims = 3, firstaxis = 1, lastaxis = 2, catype = "CA", 
%    plottype = "classic", biptype = "row", scaleplot = 3, pos = 2, 
%    M = min(nrow(Xtable), ncol(Xtable)) - 1) 
%{
%    if (!any(catype == c("CA", "SOCA", "DOCA", "NSCA", "SONSCA", 
%        "DONSCA"))) 
%        stop(paste("Must be CA, DOCA, SOCA, NSCA, SONSCA or DONSCA"))
%    X <- as.matrix(Xtable)
%    rowlabels <- rownames(Xtable)
%    collabels <- colnames(Xtable)
%    rows <- dim(X)[1]
%    cols <- dim(X)[2]
%    n <- sum(X)
%    if (is.null(mj)) {
%        mj <- c(1:cols)
%    }
%    else mj <- c(mj)
%    if (is.null(mi)) {
%        mi <- c(1:rows)
%    }
%    else mi <- c(mi)
%    maxaxes <- min(rows, cols) - 1
%    S <- switch(catype, CA = cabasic(X), SOCA = socabasic(X, 
%        mj), DOCA = docabasic(X, mi, mj), NSCA = nscabasic(X), 
%        SONSCA = sonscabasic(X, mj), DONSCA = donscabasic(X, 
%            mi, mj))
%    if (catype == "CA") {
%        Fmat <- S@RX \%*\% S@Rweights \%*\% S@Raxes
%        Gmat <- S@CX \%*\% S@Cweights \%*\% S@Caxes
%        dmum1 <- diag((S@mu + (S@mu == 0)) * (1 - (S@mu == 0)))
%        Fbi <- S@Cweights \%*\% S@Caxes
%        Gbi <- S@Rweights \%*\% S@Raxes
%        pcc <- t(S@C)
%        inertia <- S@mu * S@mu
%        comps <- diag(inertia)
%        Trend <- (Fmat[, firstaxis:lastaxis] \%*\% t(Gbi[, firstaxis:lastaxis]))
%    }
%    if (catype == "DOCA") {
%        pcc <- S@Rprofile
%        Gbi <- S@Raxes
%        Fbi <- S@Caxes
%        Gmat <- S@CX \%*\% S@Caxes
%        Fmat <- S@RX \%*\% S@Raxes
%        if ((S@Z[1, 1] < 0) || (S@Z[1, 2] < 0)) {
%            Gmat <- (-1) * Gmat
%            Fmat <- (-1) * Fmat
%        }
%        reconstruction <- t(Gmat \%*\% t(S@Cweights \%*\% Fbi))
%        dimnames(reconstruction) <- dimnames(X)
%        inertia <- S@mu
%        inertia2 <- S@mu2
%        Z <- S@Z
%        comps <- compstable.exe(Z)
%        Icompnames <- c("** Row Components **", "Location", "Dispersion", 
%            "Error", "** Column Components **", "Location", "Dispersion", 
%            "Error", "** Chi-squared Statistic **")
%        Jcompnames <- c("Component Value", "P-value")
%        dimnames(Z) <- list(paste("Poly", 1:(rows - 1)), paste("Poly", 
%            1:(cols - 1)))
%        dimnames(comps) <- list(paste(Icompnames), paste(Jcompnames))
%        Trend <- (Fmat[, firstaxis:lastaxis] \%*\% t(S@Rweights \%*\% 
%            Gbi[, firstaxis:lastaxis]))
%    }
%    if (catype == "SOCA") {
%        pcc <- S@Rprofile
%        dimnames(pcc) <- dimnames(X)
%        Gmat <- S@CX \%*\% S@Caxes
%        Fmat <- S@RX \%*\% S@Rweights \%*\% S@Raxes
%        if (S@Z[1, 1] < 0) {
%            Gmat <- (-1) * Gmat
%        }
%       Gbi <- S@Raxes
%        Fbi <- S@Caxes
%        inertia <- S@mu * S@mu
%        inertia2 <- S@mu2[-1]
%        comps <- diag(inertia2)
%        comps <- compsonetable.exe(S@Z)
%        Icompnames <- c("** Column Components **", "Location", 
%            "Dispersion", "Error", "** C-Statistic **")
%        Jcompnames <- c("Component Value", "P-value")
%        dimnames(S@Z) <- list(paste("Axis", 1:nrow(S@Z)), paste("Poly", 
%            1:(cols - 1)))
%        dimnames(comps) <- list(paste(Icompnames), paste(Jcompnames))
%        Trend <- (Fmat[, firstaxis:lastaxis] \%*\% t(S@Rweights \%*\% 
%            Gbi[, firstaxis:lastaxis]))
%    }
%    if (catype == "NSCA") {
%        Fbi <- S@Caxes
%        Gbi <- S@Raxes
%        dmum1 <- diag((S@mu + (S@mu == 0)) * (1 - (S@mu == 0)))
%        pcc <- S@RX
%        dimnames(pcc) <- dimnames(X)
%        Gmat <- S@Raxes \%*\% dmum1
%        Fmat <- S@Caxes \%*\% dmum1
%        inertia <- S@mu * S@mu
%        comps <- diag(inertia)
%        Trend <- (Fmat[, firstaxis:lastaxis] \%*\% t(S@Rweights \%*\% 
%            Gbi[, firstaxis:lastaxis]))
%    }
%    if (catype == "DONSCA") {
%        Fbi <- S@Caxes
%        Gbi <- S@Raxes
%        pcc <- S@RX
%        dimnames(pcc) <- dimnames(X)
%        Gmat <- S@CX \%*\% S@Cweights \%*\% S@Caxes
%        Fmat <- S@RX \%*\% S@Rweights \%*\% S@Raxes
%        if (S@Z[1, 1] < 0) {
%            Gmat <- (-1) * Gmat
%        }
%        inertia <- S@mu
%        inertia2 <- S@mu2
%        Z <- sqrt((n - 1) * (rows - 1)) * S@Z
%        comps <- compstable.exe(Z)
%        Icompnames <- c("** Row Components **", "Location", "Dispersion", 
%            "Error", "** Column Components **", "Location", "Dispersion", 
%            "Error", "** Chi-squared Statistic **")
%        Jcompnames <- c("Component Value", "P-value")
%        dimnames(Z) <- list(paste("Poly", 1:(rows - 1)), paste("Poly", 
%            1:(cols - 1)))
%        dimnames(comps) <- list(paste(Icompnames), paste(Jcompnames))
%        Trend <- (Fmat[, firstaxis:lastaxis] \%*\% t(Gbi[, firstaxis:lastaxis]))
%    }
%    if (catype == "SONSCA") {
%        pcc <- S@Rprofile
%        dimnames(pcc) <- dimnames(X)
%        Gmat <- S@CX \%*\% S@Cweights \%*\% S@Caxes
%        Fmat <- S@RX \%*\% (S@Rweights) \%*\% S@Raxes
%        if (S@Z[1, 1] < 0) {
%            Gmat <- (-1) * Gmat
%        }
%        Gbi <- S@Raxes
%        Fbi <- S@Caxes
%        inertia <- S@mu
%        inertia2 <- S@mu2
%        Z <- sqrt((n - 1) * (rows - 1)) * S@Z
%        comps <- compsonetable.exe(S@Z)
%        Icompnames <- c("** Column Components **", "Location", 
%            "Dispersion", "Error", "** C-Statistic **")
%        Jcompnames <- c("Component Value", "P-value")
%        dimnames(S@Z) <- list(paste("Axis", 1:nrow(S@Z)), paste("Poly", 
%            1:(cols - 1)))
%        dimnames(comps) <- list(paste(Icompnames), paste(Jcompnames))
%        Trend <- t(Gmat[, firstaxis:lastaxis] \%*\% t(Fbi[, firstaxis:lastaxis]))
%    }
%    dmum2 <- diag(1/(inertia + (S@mu == 0)) * (1 - (S@mu == 0)))
%    inertiasum <- sum(inertia)
%    inertiapc <- 100 * inertia/inertiasum
%    cuminertiapc <- cumsum(inertiapc)
%    inertiapc <- round(100 * inertiapc)/100
%    cuminertiapc <- round(100 * cuminertiapc)/100
%    inertias <- round(cbind(inertia, inertiapc, cuminertiapc), 
%        dig = 3)
%    if ((catype == "SOCA") | (catype == "SONSCA") | (catype == 
%        "DOCA") | (catype == "DONSCA")) {
%        inertiasum2 <- sum(inertia2)
%        inertiapc2 <- 100 * inertia2/inertiasum2
%        cuminertiapc2 <- cumsum(inertiapc2)
%        inertiapc2 <- round(100 * inertiapc2)/100
%        cuminertiapc2 <- round(100 * cuminertiapc2)/100
%        inertias2 <- cbind(inertia2, inertiapc2, cuminertiapc2)
%    }
%    else inertias2 <- inertias
%    Xstd <- X/sum(X)
%    if ((catype == "CA") | (catype == "SOCA") | (catype == "DOCA")) {
%        dr <- diag(rowSums(Xstd))
%    }
%    else {
%        uni <- rep(1, rows)
%        dr <- diag(uni)
%    }
%    dc <- diag(colSums(Xstd))
%    dimnames(Trend) <- list(rowlabels, collabels)
%    cacorporateplus <- new("cacorporateplus", S = S, DataMatrix = X, 
%        rows = rows, cols = cols, rowlabels = rowlabels, collabels = collabels, 
%        Rprinccoord = Fmat, Cprinccoord = Gmat, Rstdcoord = Fbi, 
%        Cstdcoord = Gbi, inertiasum = inertiasum, inertias = inertias, 
%        inertias2 = inertias2, comps = comps, printdims = printdims, 
%        maxaxes = maxaxes, catype = catype, mj = mj, mi = mi, 
%        pcc = pcc, Jmass = dc, Imass = dr, Trend = Trend)
%    print(cacorporateplus)
%    plot(cacorporateplus, cex = cex, firstaxis = firstaxis, lastaxis = lastaxis, 
%        prop = prop, M = M, catype = catype, biptype = biptype, 
%        plottype = plottype, scaleplot = scaleplot, pos = pos)
%  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{multivariate}
\keyword{multivariate}% __ONLY ONE__ keyword per line
