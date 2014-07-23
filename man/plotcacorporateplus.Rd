\name{plotcacorporateplus}
\alias{plotcacorporateplus}
\title{
Main plot function 
}
\description{
This function allows the analyst to produce the suitable graphical displays with respect to the six variants of correspondence analysis. 
The main plot function called from the main function \code{CAvariants}. 
In particular when \code{plottype = "classic"}, it produces classical graphical displays for \code{catype = "CA"} and \code{catype = "NSCA"}, 
where the row and column variables are graphically depicted in principal coordinates.\cr
When we set \code{plottype = "biplot"}, it produces biplot graphical displays, or polynomial biplots in case of ordered analysis, 
it can be specified when the input parameter \code{catype} of the main code \code{CAvariants} is equal to \code{catype = "SOCA"}, \code{catype = "DOCA"},
\code{catype = "SONSCA"} and \code{catype = "DONSCA"}. 
 }
\usage{
plotcacorporateplus(x, firstaxis = 1, lastaxis = 2, cex = 1, 
prop = 1, inert, inertsum, Xdata, catype = "CA", biptype = "rows", 
plottype = "classic", M = (min(nrow(x), ncol(x)) - 1), 
scaleplot = 3, posleg = "topleft", pos = 2, ell = TRUE)
}
\arguments{
  \item{x}{
Represents the set of the output parameters of the main function \code{CAvariants} of the \code{R} object class \code{cacorporateplus}.
}
  \item{firstaxis}{
The dimension reflected along  the horizontal axis.
}
  \item{lastaxis}{
The dimension reflected along the vertical axis.
}
  \item{cex}{
The parameter that specifies the size of character labels of points in graphical displays. By default, it is equal to 1.
}
  \item{prop}{
The scaling parameter for specifying the limits of the plotting area. By default, it is equal to 1. 
}
  \item{inert}{
 The percentage of explained inertia for each axis.
}
  \item{inertsum}{
 The total amount of point inertia.
}
  \item{Xdata}{
The two-way contingency table.
}
  \item{catype}{
 The input parameter specifying what variant of correspondence analysis is requested.
}
  \item{plottype}{
 The input parameter specifying what kind of graphical display is requested. By default, it is equal 
to \code{plottype = "classic"}.
}
\item{biptype}{
The input parameter specifying what kind of biplot is requested. By default, it is equal to \code{biptype = "row"}
(row isometric biplot).
}
\item{M}{
The rank of the two-way contingency table.
}
\item{scaleplot}{
The scaling parameter for scaling biplot coordinates. By default, \code{scaleplot = 1}.
}
\item{posleg}{
The parameter for specifying the position of the legend in the graphical function \code{trendplot}. 
By default, \code{posleg= "topleft"}.}
 \item{pos}{
The parameter for specifying the position of  point symbols in the graphical displays \cr
By default, \code{pos=2}.}
\item{ell}{ The flag parameter to allow the building of algebraic ellipses. By default, \code{ell = T}.}
}
\details{
It produces classical and biplot graphical displays. Further when \code{catype} is equal to \cr
\code{"DOCA", "SOCA", "DONSCA"} or \code{"SONSCA"}, 
the trend of row and column variables after the reconstruction of column profiles by polynomials is portrayed.\cr
For classical biplot displays, it superimposes on it algebraic ellipses of confidence. It uses the secondary plot functions \code{caellipse} or
\code{nsca.ellipse}, depending on the input parameter \code{catype}. }
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
Rosaria Lombardo and Eric J Beh
}
\note{
For classical graphical displays, both sets of coordinates are defined using principal coordinates (see Greenacre's terminology).
In biplot graphical displays, one set of coordinates is standard and the other is principal. 
When the analysis is ordered, one set of coordinates consists of standard polynomial coordinates and the other one is of principal polynomial coordinates.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%\examples{
%plot.cacorporateplus(x)
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%## The function is currently defined as
%function (x, firstaxis = 1, lastaxis = 2, cex = cex, 
%    picsize = c(-1, 1), prop = 1, inert, inertsum, Xdata, catype = "CA", 
%    biptype = "rows", plottype = "classic", M = (min(nrow(x), 
%        ncol(x)) - 1), scaleplot = 3, pos = 2) 
%{
%    if ((firstaxis < 1) | (firstaxis > x@maxaxes - 1)) 
%        stop(paste("incorrect first axis =", firstaxis, "\n\n"))
%    if (lastaxis > x@maxaxes) 
%        stop(paste("incorrect last axis =", lastaxis, "\n\n"))
%    if (firstaxis >= lastaxis) 
%        stop(paste("last axis must be greater than first axis\n\n"))
%    if (dim(array(picsize)) != 2) 
%        stop(paste("need lower and upper bound picsize =", picsize, 
%            "\n\n"))
%    if (picsize[1] >= picsize[2]) 
%        stop(paste("incorrect axis scale picsize =", picsize[1], 
%            picsize[2], "\n\n"))
%    tworowS <- rowSums(x@DataMatrix > 0) == 2
%    twocolS <- colSums(x@DataMatrix > 0) == 2
%    if (is.null(groupfile)) {
%        rowgroup <- list(1:x@rows, rep(1, x@rows))
%        rowgrlab <- list(1, "", "*", "blue", "T")
%        colgroup <- list(1:x@cols, rep(1, x@cols))
%        colgrlab <- list(1, "", "+", "red", "T")
%    }
%    else {
%        rowgroup <- scan(file = groupfile, what = list("integer", 
%            "integer"), nmax = x@rows)
%        nrowgroups <- as.integer(max(rowgroup[[2]]))
%        rowgrlab <- scan(file = groupfile, what = list("integer", 
%            "character", "character", "character", "logical"), 
%            nmax = nrowgroups, skip = x@rows + 1)
%        colgroup <- scan(file = groupfile, what = list("integer", 
%            "integer"), nmax = x@cols, skip = x@rows + nrowgroups + 
%            2)
%        ncolgroups <- as.integer(max(colgroup[[2]]))
%        colgrlab <- scan(file = groupfile, what = list("integer", 
%            "character", "character", "character", "logical"), 
%            nmax = ncolgroups, skip = x@rows + nrowgroups + x@cols + 
%                3)
%    }
%    if ((plottype == "Classical") | (plottype == "classical") | 
%        (plottype == "classic") | (plottype == "c")) {
%        nthings <- x@cols
%        nvars <- x@rows
%        cord1 <- x@Cprinccoord
%        cord2 <- x@Rprinccoord
%        inertiapc = x@inertias[, 2]
%        if ((x@catype == "DONSCA") | (x@catype == "DOCA") | (x@catype == 
%            "SOCA") | (x@catype == "SONSCA")) {
%            cord1 <- x@Cstdcoord/scaleplot
%            cord2 <- x@Rprinccoord * scaleplot
%            inertiapc = x@inertias2[, 2]
%        }
%        dimnames(cord1)[1] <- dimnames(x@DataMatrix)[2]
%        dimnames(cord2)[1] <- dimnames(x@DataMatrix)[1]
%        thinggroup <- colgroup
%        thinggrlab <- rowgrlab
%        vargroup <- rowgroup
%        vargrlab <- colgrlab
%        thinglabels <- x@collabels
%        varlabels <- x@rowlabels
%        main = "Classical plot"
%    }
%    if ((plottype == "Biplot") | (plottype == "biplot") | (plottype == 
%        "bip") | (plottype == "b")) {
%        if ((biptype == "rows") | (biptype == "Rows") | (biptype == 
%            "row") | (biptype == "r")) {
%            plottype <- "biplot"
%            cord2 <- x@Rprinccoord * scaleplot
%            cord1 <- x@Cstdcoord/scaleplot
%            nthings <- x@cols
%            nvars <- x@rows
%            thinglabels <- x@collabels
%            varlabels <- x@rowlabels
%            thinggroup <- colgroup
%            thinggrlab <- colgrlab
%            vargroup <- rowgroup
%            vargrlab <- rowgrlab
%            main <- "Row Isometric Biplot"
%            inertiapc = x@inertias[, 2]
%            if ((x@catype == "DONSCA") | (x@catype == "DOCA") | 
%                (x@catype == "SOCA") | (x@catype == "SONSCA")) {
%                inertiapc = x@inertias2[, 2]
%            }
%            dimnames(cord1)[1] <- dimnames(x@DataMatrix)[2]
%            dimnames(cord2)[1] <- dimnames(x@DataMatrix)[1]
%        }
%        else {
%            plottype <- "biplot"
%            cord1 <- x@Cprinccoord * scaleplot
%            cord2 <- x@Rstdcoord/scaleplot
%            nthings <- x@cols
%            nvars <- x@rows
%            thinggroup <- colgroup
%            thinggrlab <- colgrlab
%            vargroup <- colgroup
%            vargrlab <- rowgrlab
%            thinglabels <- x@collabels
%            varlabels <- x@rowlabels
%            main <- "Column Isometric Biplot"
%            inertiapc = x@inertias[, 2]
%            if ((x@catype == "DONSCA") | (x@catype == "DOCA") | 
%                (x@catype == "SOCA") | (x@catype == "SONSCA")) {
%                inertiapc = x@inertias2[, 2]
%            }
%            dimnames(cord1)[1] <- dimnames(x@DataMatrix)[2]
%            dimnames(cord2)[1] <- dimnames(x@DataMatrix)[1]
%        }
%    }
%    if ((catype == "DOCA") | (catype == "SOCA") | (catype == 
%        "SONSCA") | (catype == "DONSCA")) {
%        cat("\n Looking at the Trends of rows and columns\n")
%        trendplot(x@mj, (x@Trend), position = "bottomleft", main = "Reconstructed rows of the Centred column profile", 
%            xlab = "ordered scores", prop = prop)
%        windows()
%        trendplot(x@mi, t(x@Trend), position = "bottomleft", 
%            main = "Reconstrucetd columns of the Centred column profile", 
%            xlab = "ordered scores", prop = prop)
%    }
%    picsize1 <- c(range(cord1[, c(firstaxis, lastaxis)], cord2[, 
%        c(firstaxis, lastaxis)])/prop)
%    picsizeR <- c(range(cord1[, c(firstaxis, lastaxis)], cord1[, 
%        c(firstaxis, lastaxis)])/prop)
%    picsize2 <- c(range(cord2[, c(firstaxis, lastaxis)], cord2[, 
%        c(firstaxis, lastaxis)])/prop)
%    if (picsize1[1] >= picsize1[2]) 
%        stop(paste("incorrect axis scale picsize =", picsize1[1], 
%            picsize1[2], "\n\n"))
%    if ((x@catype == "DONSCA") || (x@catype == "DOCA") || (x@catype == 
%        "SOCA") || (x@catype == "SONSCA")) {
%        plotone(firstaxis, lastaxis, plottype = plottype, things = x@catype, 
%            nthings, nvars, cord1, cord2, inertiapc = inertiapc, 
%            thinggroup, thinggrlab, vargroup, vargrlab, thinglabels, 
%            varlabels, picsize = picsize1, cex = cex, type = "b", 
%            catype = "DOCA", pos = pos)
%    }
%    else {
%        plotone(firstaxis, lastaxis, plottype = plottype, things = x@catype, 
%            nthings, nvars, cord1, cord2, inertiapc = inertiapc, 
%            thinggroup, thinggrlab, vargroup, vargrlab, thinglabels, 
%            varlabels, picsize = picsize1, cex = cex, type = "p", 
%            catype = "CA", pos = pos)
%    }
%    windows()
%    M <- (min(nrow(x@DataMatrix), ncol(x@DataMatrix)) - 1)
%    switch(x@catype, CA = caellipse(x@DataMatrix, M = M, prop = prop, 
%        Imass = x@Imass, Jmass = x@Jmass, a = x@S@Caxes, b = x@S@Raxes, 
%        g = cord1, f = cord2, dmu = diag(x@inertias[, 1]), inertiapc = x@inertias[, 
%            2], pos = pos), SOCA = caellipse(x@DataMatrix, M = M, 
%        prop = prop, Imass = x@Imass, Jmass = x@Jmass, a = x@S@Caxes, 
%        b = x@S@Raxes, g = cord1, f = cord2, dmu = diag(x@inertias2[, 
%            1]), inertiapc = x@inertias2[, 2], pos = pos), DOCA = caellipse(x@DataMatrix, 
%        M = (M - 1), prop = prop, Imass = x@Imass, Jmass = x@Jmass, 
%        a = x@S@Caxes, b = x@S@Raxes, f = cord2, g = cord1, dmu = diag(x@inertias2[, 
%            1]), inertiapc = x@inertias2[, 2], pos = pos), NSCA = nsca.ellipse(x@DataMatrix, 
%        M = M, prop = prop, Imass = x@Imass, Jmass = x@Jmass, 
%        a = x@S@Caxes, b = x@S@Raxes, f = cord2, g = cord1, dmu = diag(x@inertias[, 
%            1]), inertiapc = x@inertias[, 2], plottype = plottype, 
%        catype = "NSCA", pos = pos), SONSCA = nsca.ellipse(x@DataMatrix, 
%        M = M, prop = prop, Imass = x@Imass, Jmass = x@Jmass, 
%        a = x@S@Caxes, b = x@S@Raxes, f = cord2, g = cord1, dmu = diag(x@inertias2[, 
%            1]), inertiapc = x@inertias2[, 2], plottype = plottype, 
%        catype = "SONSCA", pos = pos), DONSCA = nsca.ellipse(x@DataMatrix, 
%        M = (M - 1), prop = prop, Imass = x@Imass, Jmass = x@Jmass, 
%        a = x@S@Caxes, b = x@S@Raxes, f = cord2, g = cord1, dmu = diag(x@inertias2[, 
%            1]), inertiapc = x@inertias2[, 2]), plottype = plottype, 
%        catype = "DONSCA", pos = pos)
%  }
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
\keyword{multivariate}% __ONLY ONE__ keyword per line
