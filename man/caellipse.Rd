\name{caellipse}
\alias{caellipse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Algebraic elliptical confidence regions 
}
\description{
It produces elliptical confidence regions when symmetrical or ordered symmetrical correspondence analysis is performed. \cr
This function allows the analyst to superimpose confidence ellipses into a graphical display when the input parameter \code{catype}
of the main code \code{CAvariants} is equal to \code{CA, SOCA} or \code{DOCA}. \cr
It is called internally  from the main plot function \code{plot.CAvariants}. \cr
It uses the function \code{ellipse}.

}
\usage{
caellipse(N, a1 = 1, a2 = 2, alpha = 0.05, cols = c(2, 4), 
M = 2, cex = 0.8, cex.lab = 0.5,  
mar = c(5, 4, 4, 2) + 0.1, prop = 0.8, Inames, Jnames, Imass, 
Jmass, a, b, g, f, dmu, inertiapc, pos = 2, plottype="biplot", 
biptype = "row", arrow = TRUE, length = 0.01)
}
\arguments{
\item{N}{
The two-way contingency table.
}
\item{a1}{
The dimension reflected along the horizontal axis.
}
  \item{a2}{
The dimension reflected along the vertical axis.
}
  \item{alpha}{
The confidence level of the elliptical regions.
}
  \item{cols}{
The graphical parameter for setting the colours of points in the graphical displays.
}
  \item{M}{
The number of  axes used when drawing ellypse. By default, it is equal to \code{M= 2}.
Its maximum value is equal to the rank of the data matrix.}
  \item{cex}{
The parameter for setting the size of character labels of points in graphical displays.
}
  \item{cex.lab}{
The parameter for setting the size of character labels of axes in graphical displays.
}
  \item{mar}{
The parameter for setting the size of the plotting area. 
}
  \item{prop}{
The scaling parameter for specifying the limits of the plotting area. By default it is equal to 1. 
}
  \item{Inames}{
The labels of the row variable.}
  \item{Jnames}{
The labels of the column variable.
}
  \item{Imass}{
The weight matrix of the row variable.
}
  \item{Jmass}{
The weight matrix of the column variable.
}
  \item{a}{
The row principal or polynomial axes.
}
  \item{b}{
The column principal or polynomial axes.
}
  \item{f}{
The row coordinates.
}
  \item{g}{
The column coordinates.
}
  \item{dmu}{
The squared singular values or inertia of each axis.
}
  \item{inertiapc}{
The percentage of explained inertia.
}
  \item{pos}{
The parameter for specifying the position of  point symbols in the graphical displays. \cr
By default, \code{pos=2}.}
  \item{plottype}{
The parameter for specifying the type of plot. By default, \code{plottype="classic"}. 
}
 \item{biptype}{
The parameter for specifying the type of biplot. By default, \code{biptype="row"}. 
}
\item{arrow}{
The parameter for diplaying arrows in biplot in correspondence of points in principal coordinates. By default, \code{arrow=TRUE}. 
}
\item{length}{
The parameter for setting the length of the arrow end in biplot. By default, \code{length=0.01}. 
}

}
\details{
The output values of this function
}
\value{
  \item{eccentricity}{Value of ellipse eccentricity.}
  \item{row.summ }{Value of ellipse semi-axes for each row points.}
  \item{col.summ }{Value of ellipse semi-axes for each column points.}
}
\references{
Beh EJ and Lombardo R 2014 Correspondence analysis, Theory, Practice and New Strategies. John Wiley & Sons.
}
\author{
Rosaria Lombardo and Eric J Beh
}
\note{
This function is called from the main plot function \code{plot.CAvariants} and is executed when \code{catype="CA"}, \code{catype="SOCA"} or \code{"DOCA"}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%\examples{
%caellipse(asbestos, a1 = 1, a2 = 2, alpha = 0.05, cols = c(2, 4), M = min(nrow(N), ncol(N)) - 1, cex = 0.8, cex.lab = 0.5, 
%mar = c(5, 4, 4, 2) + 0.1, prop = 0.8, Inames, Jnames, Imass, Jmass, a, b, f, g, dmu, inertiapc, posl = 2)
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%## The function is currently defined as
%function (N, a1 = 1, a2 = 2, alpha = 0.05, cols = c(2, 4), M = min(nrow(N), 
%    ncol(N)) - 1, cex = 0.8, cex.lab = 0.5, mar = c(5, 4, 4, 
%    2) + 0.1, n1 = 1, n2 = 2, prop = 0.8, Inames, Jnames, Imass, 
%    Jmass, a, b, f, g, dmu, inertiapc, pos = 2) 
%{
%    I <- nrow(f)
%    J <- nrow(g)
%    n <- sum(N)
%    rowgroup <- list(1:I, rep(1, I))
%    rowgrlab <- list(1, "", "*", "blue", "T")
%    colgroup <- list(1:J, rep(1, J))
%    colgrlab <- list(1, "", "+", "red", "T")
%    Inames <- dimnames(f)[1]
%    Jnames <- dimnames(g)[1]
%    dIh <- solve(Imass^0.5)
%    dJh <- solve(Jmass^0.5)
%    t.inertia <- sum(dmu)/n
%    dmu <- sqrt(dmu/n)
%    a <- dIh \%*\% a
%    b <- dJh \%*\% b
%    chisq.val <- qchisq(1 - alpha, df = (I - 1) * (J - 1))
%    hlax1.row <- vector(mode = "numeric", length = I)
%    hlax2.row <- vector(mode = "numeric", length = I)
%    hlax1.col <- vector(mode = "numeric", length = J)
%    hlax2.col <- vector(mode = "numeric", length = J)
%    if (M > 2) {
%        for (i in 1:I) {
%            hlax1.row[i] <- dmu[1, 1] * sqrt(abs((chisq.val/(n * 
%                t.inertia)) * (1/Imass[i, i] - sum(a[i, 3:M]^2))))
%            hlax2.row[i] <- dmu[2, 2] * sqrt(abs((chisq.val/(n * 
%                t.inertia)) * (1/Imass[i, i] - sum(a[i, 3:M]^2))))
%        }
%        for (j in 1:J) {
%            hlax1.col[j] <- dmu[1, 1] * sqrt(abs((chisq.val/(n * 
%                t.inertia)) * (1/Jmass[j, j] - sum(b[j, 3:M]^2))))
%            hlax2.col[j] <- dmu[2, 2] * sqrt(abs((chisq.val/(n * 
%                t.inertia)) * (1/Jmass[j, j] - sum(b[j, 3:M]^2))))
%        }
%    }
%    else {
%        for (i in 1:I) {
%            hlax1.row[i] <- dmu[1, 1] * sqrt(abs((chisq.val/(n * 
%                t.inertia)) * (1/(Imass)[i, i])))
%            hlax2.row[i] <- dmu[2, 2] * sqrt(abs((chisq.val/(n * 
%                t.inertia)) * (1/(Imass)[i, i])))
%        }
%        for (j in 1:J) {
%            hlax1.col[j] <- dmu[1, 1] * sqrt(abs((chisq.val/(n * 
%                t.inertia)) * (1/(Jmass)[j, j])))
%            hlax2.col[j] <- dmu[2, 2] * sqrt(abs((chisq.val/(n * 
%                t.inertia)) * (1/(Jmass)[j, j])))
%        }
%    }
%    plot.new()
%    par(pty = "s", mar = mar)
%    plot(0, 0, xlim = range(f, g)/prop, ylim = range(f, g)/prop, 
%        xlab = paste("Axis ", a1, "    ", inertiapc[a1], "\%", 
%            sep = ""), ylab = paste("Axis ", a2, "    ", inertiapc[a2], 
%            "\%", sep = ""), asp = 1, pch = rowgrlab[[3]][as.integer(rowgroup[[2]])], 
%        col = rowgrlab[[4]][as.integer(rowgroup[[2]])], cex = cex)
%    text(f[, a1], f[, a2], labels = Inames[[1]], pos = pos, col = cols[1], 
%        cex = cex)
%    points(f[, a1], f[, a2], pch = "*", col = cols[1])
%    text(g[, a1], g[, a2], labels = Jnames[[1]], pos = 2, col = cols[2], 
%        cex = cex)
%    points(g[, a1], g[, a2], pch = "+", col = cols[2])
%    abline(h = 0, v = 0)
%    title(main = paste(100 * (1 - alpha), "\% Confidence Ellipses"))
%    for (i in 1:I) {
%        ellipse(hlax1.row[i], hlax2.row[i], xc = f[i, 1], yc = f[i, 
%            2], col = cols[1])
%    }
%    for (j in 1:J) {
%        ellipse(hlax1.col[j], hlax2.col[j], , xc = g[j, 1], yc = g[j, 
%            2], col = cols[2])
%    }
%    eccentricity <- sqrt(1 - (dmu[2, 2]/dmu[1, 1])^2)
%    area.row <- vector(mode = "numeric", length = I)
%    area.col <- vector(mode = "numeric", length = J)
%    for (i in 1:I) {
%        area.row[i] <- 3.14159 * hlax1.row[i] * hlax2.col[i]
%    }
%    for (j in 1:J) {
%        area.col[j] <- 3.14159 * hlax1.col[j] * hlax2.col[j]
%    }
%    pvalrow <- vector(mode = "numeric", length = I)
%    pvalcol <- vector(mode = "numeric", length = J)
%    for (i in 1:I) {
%        if (M > 2) {
%            pvalrow[i] <- 1 - pchisq(n * t.inertia * ((1 - sum(a[i, 
%                3:M]^2))^(-1)) * (f[i, 1]^2/dmu[1, 1]^2 + f[i, 
%                2]^2/dmu[2, 2]^2), df = (I - 1) * (J - 1))
%        }
%        else {
%            pvalrow[i] <- 1 - pchisq(n * t.inertia * ((1)^(-1)) * 
%                (f[i, 1]^2/dmu[1, 1]^2 + f[i, 2]^2/dmu[2, 2]^2), 
%                df = (I - 1) * (J - 1))
%        }
%    }
%    for (j in 1:J) {
%        if (M > 2) {
%            pvalcol[j] <- 1 - pchisq(n * t.inertia * ((1/Jmass[j, 
%                j] - sum(b[j, 3:M]^2))^(-1)) * (g[j, 1]^2/dmu[1, 
%                1]^2 + g[j, 2]^2/dmu[2, 2]^2), df = (I - 1) * 
%                (J - 1))
%        }
%        else {
%            pvalcol[j] <- 1 - pchisq(n * t.inertia * ((1/Jmass[j, 
%                j])^(-1)) * (g[j, 1]^2/dmu[1, 1]^2 + g[j, 2]^2/dmu[2, 
%                2]^2), df = (I - 1) * (J - 1))
%        }
%    }
%    summ.name <- c("HL Axis 1", "HL Axis 2", "Area", "P-value")
%    row.summ <- cbind(hlax1.row, hlax2.row, area.row, pvalrow)
%    dimnames(row.summ) <- list(paste(Inames[[1]]), paste(summ.name))
%    col.summ <- cbind(hlax1.col, hlax2.col, area.col, pvalcol)
%    dimnames(col.summ) <- list(paste(Jnames[[1]]), paste(summ.name))
%    invisible(list(eccentricity = eccentricity, row.summ = row.summ, 
%        col.summ = col.summ))
%  }
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
