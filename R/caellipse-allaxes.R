caellipse <-
function (Xtable, a1=1,a2=2,alpha = 0.05, cols = c(2, 4), M = 2, cex =0.8, cex.lab = 0.8, mar = c(5, 4, 4, 2) + 0.1, 
 prop=0.8,Imass,Jmass,a,b,g,fr,dmu,inertiapc,plottype="biplot",biptype="row",
pos=2,arrow=TRUE,length=0,graphy=TRUE,ell=TRUE) 
{
    I <- ncol(Imass)
    J <- ncol(Jmass)
n<-sum(Xtable)
MI<-ncol(a)
MJ<-ncol(b)
rowgroup <- list(1:I,rep(1,I))
rowgrlab <- list(1,"","*","blue","T")
colgroup <- list(1:J,rep(1,J))
colgrlab <- list(1,"","+","red","T")
    Inames <- dimnames(fr)[1]
    Jnames <- dimnames(g)[1]
    dIh <- solve(Imass^0.5)
    dJh <- solve(Jmass^0.5)
    t.inertia <- sum(dmu)
dmu<-sqrt(dmu)
chisq.val <- qchisq(1 - alpha, df = (MI - 1) * (MJ - 1))
    hlax1.row <- vector(mode = "numeric", length = I)
    hlax2.row <- vector(mode = "numeric", length = I)
    hlax1.col <- vector(mode = "numeric", length = J)
    hlax2.col <- vector(mode = "numeric", length = J)
#browser()
    if (M > 2) {
        for (i in 1:I) {
            hlax1.row[i] <- sqrt(dmu[1,1]) * sqrt(abs((chisq.val/(t.inertia*n)) * 
                (1/Imass[i,i] - sum(a[i, 3:MI ]^2))))
            hlax2.row[i] <- sqrt(dmu[2,2]) * sqrt(abs((chisq.val/(t.inertia*n)) * 
                (1/Imass[i,i] - sum(a[i, 3:MI]^2))))
        }
        for (j in 1:J) {
            hlax1.col[j] <- sqrt(dmu[1,1] )* sqrt(abs((chisq.val/(t.inertia*n)) * 
                (1/Jmass[j,j] - sum(b[j, 3:MJ]^2))))
            hlax2.col[j] <- sqrt(dmu[2,2]) * sqrt(abs((chisq.val/(t.inertia*n)) * 
                (1/Jmass[j,j] - sum(b[j, 3:MJ]^2))))
        }
    }
    else {
        for (i in 1:I) {
            hlax1.row[i] <- sqrt(dmu[1,1]) * sqrt(abs((chisq.val/(t.inertia*n)) * 
                ((1/Imass)[i,i])))
            hlax2.row[i] <- sqrt(dmu[2,2]) * sqrt(abs((chisq.val/(t.inertia*n)) * 
                ((1/Imass)[i,i])))
        }
        for (j in 1:J) {
            hlax1.col[j] <- sqrt(dmu[1,1]) * sqrt(abs((chisq.val/(t.inertia*n)) * 
                ((1/Jmass)[j,j])))
            hlax2.col[j] <- sqrt(dmu[2,2]) * sqrt(abs((chisq.val/(t.inertia*n)) * 
                ((1/Jmass)[j,j])))
        }
    }
if (graphy==TRUE){
picsize<-c(range(fr[,c(a1,a2)], g[,c(a1,a2)])/prop)
    plot.new()
    par(pty = "s", mar = mar)
grat <- cbind(fr[,a1]/picsize[1],fr[,a1]/picsize[2],fr[,a2]/picsize[1],fr[,a2]/picsize[2],0.95)/0.95
 cl <- 1.05/apply(grat,1,max)
grat2 <- cbind(g[,a1]/picsize[1],g[,a1]/picsize[2],g[,a2]/picsize[1],g[,a2]/picsize[2],0.95)/0.95
  cl2 <- 1.05/apply(grat2,1,max)
#browser()
#plot(fr[, a1], fr[, a2], xlim=picsize, ylim=picsize, 
plot(0,0,pch=" ", xlim=picsize, ylim=picsize, 
 xlab=paste("Axis ", a1, "    ", inertiapc[a1], "%", sep=""), 
 ylab=paste("Axis ", a2, "    ", inertiapc[a2], "%", sep=""),  asp=1,col=rowgrlab[[4]][as.integer(rowgroup[[2]])],cex=cex,cex.lab=cex.lab)
abline(h = 0, v = 0) 
text(cl*fr[, a1], cl*fr[, a2], labels = Inames[[1]], pos= pos, col = cols[1], cex = cex)
points(cl*fr[, a1],cl* fr[, a2], pch="*",col=cols[1],cex=cex)
    text(cl2*g[, a1],cl2* g[, a2], labels = Jnames[[1]], pos= pos, col = cols[2], 
        cex = cex)
points(cl2*g[, a1], cl2*g[, a2], pch="+",col=cols[2],cex=cex)
if (ell){    
   title(main = paste(100 * (1 - alpha), "% Confidence Ellipses"))}
else {title(main =  "Row Biplot")}
if (plottype=="biplot") {
if ((biptype=="column")|(biptype=="col")|(biptype=="c")) {
#----------------------------------------------arrow on column principal coords
nv <- rep(0, nrow(g))
#vec <- g[, c(a1,a2)]
if(arrow) {
arrows(nv, nv, g[,a1], g[, a2], length = length)
}
 for (j in 1:J) {
        ellipse(hlax1.col[j], hlax2.col[j], xc = g[j, a1], yc = g[j, 
            a2], col = cols[2])
    }
} #end col
if ((biptype=="row")|(biptype=="r")|(biptype=="rows")){
#----------------------------------------------arrow on row principal coords
nv <- rep(0, nrow(fr))
if(arrow) {
arrows(nv, nv, fr[, a1], fr[,a2], length = length)
}
if (ell){
 for (i in 1:I) {
        ellipse(hlax1.row[i], hlax2.row[i], xc = fr[i, a1], yc = fr[i, 
            a2], col = cols[1])
    }
} }#end row end ell
}#end biplot
else{ # when classical 
 for (j in 1:J) {
        ellipse(hlax1.col[j], hlax2.col[j], xc = g[j, a1], yc = g[j, 
            a2], col = cols[2])
    }
for (i in 1:I) {
        ellipse(hlax1.row[i], hlax2.row[i], xc = fr[i, a1], yc = fr[i, 
           a2], col = cols[1])
    }
}#end else
}#end graphy
    eccentricity <- abs(1 - (dmu[2, 2]^2/dmu[1, 1]^2))^(1/2)
    area.row <- vector(mode = "numeric", length = I)
    area.col <- vector(mode = "numeric", length = J)
    for (i in 1:I) {
        area.row[i] <- 3.14159 * hlax1.row[i] * hlax2.row[i]
    }
    for (j in 1:J) {
        area.col[j] <- 3.14159 * hlax1.col[j] * hlax2.col[j]
    }
    pvalrow <- vector(mode = "numeric", length = I)
    pvalcol <- vector(mode = "numeric", length = J)
    for (i in 1:I) {
        if (M > 2) {
            pvalrow[i] <- 1- pchisq(t.inertia *n* ((1/Imass[i,i] - sum(a[i, 
                3:MI]^2))^(-1)) * (fr[i, 1]^2/dmu[1, 1] + fr[i, 
                2]^2/dmu[2, 2]), df =(I-1)*(J-1) )
        }
        else {
            pvalrow[i] <- 1-pchisq(t.inertia*n * (1/Imass[i,i]) * 
                (fr[i, 1]^2/dmu[1, 1] + fr[i, 2]^2/dmu[2, 2]), 
                df =  (J-1) )
        }
    }
    for (j in 1:J) {
        if (M > 2) {
            pvalcol[j] <-  1-pchisq(t.inertia *n* ((1/Jmass[j,j] - 
                sum(b[j, 3:MJ]^2))^(-1)) * (g[j, 1]^2/dmu[1, 1] + 
                g[j, 2]^2/dmu[2, 2]), df = (I-1)*(J-1)  )
        }
        else {
            pvalcol[j] <- 1- pchisq(t.inertia *n* (1/Jmass[j,j]) * 
                (g[j, 1]^2/dmu[1, 1] + g[j, 2]^2/dmu[2, 2]), 
                df = (J-1) )
        }
    }
    summ.name <- c("HL Axis 1", "HL Axis 2", "Area", "P-value")
    row.summ <- cbind(hlax1.row, hlax2.row, area.row, pvalrow)
#browser()
    dimnames(row.summ) <- list(paste(Inames[[1]]), paste(summ.name))
    col.summ <- cbind(hlax1.col, hlax2.col, area.col, pvalcol)
    dimnames(col.summ) <- list(paste(Jnames[[1]]), paste(summ.name))
#browser() 
 invisible(list(eccentricity = eccentricity, row.summ = row.summ, 
        col.summ = col.summ))
}
