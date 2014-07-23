nsca.ellipse <-
function(N, a1=1,a2=2,alpha = 0.05, cols = c(2, 4), M = min(nrow(N), ncol(N))-1, cex = .8, cex.lab = 0.5, mar = c(5, 4, 4, 2) + 0.1,  
prop=.8,Inames,Jnames,Imass,Jmass,a,b,f,g,dmu,inertiapc,plottype="biplot",biptype="row",pos=2,arrow=TRUE,length=0.01) 
 {
	  I <- nrow(f)
    J <- nrow(g)
	n<-sum(N)
rowgroup <- list(1:I,rep(1,I))
rowgrlab <- list(1,"","*","blue","T")
colgroup <- list(1:J,rep(1,J))
colgrlab <- list(1,"","+","red","T")
    Inames <- dimnames(f)[1]
    Jnames <- dimnames(g)[1]
    dIh <- solve(Imass^0.5)
    dJh <- solve(Jmass^0.5)
   dmu<-sqrt(dmu)
        t.inertia <- sum(dmu^2)
    chisq.val <- (qchisq(1 - alpha, df = (I - 1) * (J - 1)))/(n - 1) * (I - 1)
    hlax1.row <- vector(mode = "numeric", length = I)
    hlax2.row <- vector(mode = "numeric", length = I)
    hlax1.col <- vector(mode = "numeric", length = J)
    hlax2.col <- vector(mode = "numeric", length = J)
if (M > 2) {
        for (i in 1:I) {
            hlax1.row[i] <- dmu[1, 1] * sqrt(abs((chisq.val/t.inertia) * 
                (1 - sum(a[i, 3:M]^2))))
            hlax2.row[i] <- dmu[2, 2] * sqrt(abs((chisq.val/t.inertia) * 
                (1 - sum(a[i, 3:M]^2))))
        }
        for (j in 1:J) {
            hlax1.col[j] <- dmu[1, 1] * sqrt(abs((chisq.val/t.inertia) * 
                (1/Jmass[j,j] - sum(b[j, 3:M]^2))))
            hlax2.col[j] <- dmu[2, 2] * sqrt(abs((chisq.val/t.inertia) * 
                (1/Jmass[j,j] - sum(b[j, 3:M]^2))))
        }
    }
    else {
        for (i in 1:I) {
            hlax1.row[i] <- dmu[1, 1] * sqrt(abs((chisq.val/t.inertia) * 
                (1)))
            hlax2.row[i] <- dmu[2, 2] * sqrt(abs((chisq.val/t.inertia) * 
                (1)))
        }
        for (j in 1:J) {
            hlax1.col[j] <- dmu[1, 1] * sqrt(abs((chisq.val/t.inertia) * 
                (1/Jmass[j,j])))
            hlax2.col[j] <- dmu[2, 2] * sqrt(abs((chisq.val/t.inertia) * 
                (1/Jmass[j,j])))
        }
    }
   #windows()
    par(pty = "s", mar = mar)
   
plot(0, 0, xlim=range(f, g)/prop, ylim=range(f, g)/prop, 
     xlab=paste("Axis ", a1, "    ", inertiapc[a1], "%", sep=""), 
     ylab=paste("Axis ", a2, "    ", inertiapc[a2], "%", sep=""),  
 asp=1, pch=rowgrlab[[3]][as.integer(rowgroup[[2]])], col=rowgrlab[[4]][as.integer(rowgroup[[2]])],cex=cex)
text(f[, 1], f[, 2], labels = Inames[[1]],pos=pos,col = cols[1], 
        cex = cex)
points(f[, a1], f[, a2], pch="+",col=cols[1])
   
 text(g[, 1], g[, 2], labels = Jnames[[1]], pos=pos, col = cols[2], 
        cex = cex)
points(g[, a1], g[, a2], pch="*",col=cols[2])

    abline(h = 0, v = 0)
    title(main = paste(100 * (1 - alpha), "% Confidence Ellipses"))
    if (plottype=="classic"){

 for (j in 1:J) {
        ellipse(hlax1.col[j], hlax2.col[j], , xc = g[j, 1], yc = g[j, 
            2], col = cols[2])
    }
 
for (i in 1:I) {
        ellipse(hlax1.row[i], hlax2.row[i], xc = f[i, 1], yc = f[i, 
            2], col = cols[1])
    }
}
else {
if ((biptype=="column")|(biptype=="col")|(biptype=="c")){
#----------------------------------------------arrow on column principal coords
nv <- rep(0, nrow(g))
	vec <- g[, c(1, 2)]

	if(arrow) {

		arrows(nv, nv, vec[, 1], vec[, 2], length = length)
	}
	
 for (j in 1:J) {
        ellipse(hlax1.col[j], hlax2.col[j], , xc = g[j, 1], yc = g[j, 
            2], col = cols[2])
    }
}
if (biptype=="row"){
#----------------------------------------------arrow on row principal coords
nv <- rep(0, nrow(f))
	vec <- f[, c(1, 2)]

	if(arrow) {

		arrows(nv, nv, vec[, 1], vec[, 2], length = length)
	}
	

#-------------------------------------------------------------
 
for (i in 1:I) {
        ellipse(hlax1.row[i], hlax2.row[i], xc = f[i, 1], yc = f[i, 
            2], col = cols[1])
    }
}
}#end else
   eccentricity <- sqrt(abs(1 - (dmu[2, 2]/dmu[1, 1])^2))
    area.row <- vector(mode = "numeric", length = I)
    area.col <- vector(mode = "numeric", length = J)
    for (i in 1:I) {
        area.row[i] <- 3.14159 * hlax1.row[i] * hlax2.col[i]
    }
    for (j in 1:J) {
        area.col[j] <- 3.14159 * hlax1.col[j] * hlax2.col[j]
    }
    pvalrow <- vector(mode = "numeric", length = I)
    pvalcol <- vector(mode = "numeric", length = J)
    for (i in 1:I) {
        if (M > 2) {
            pvalrow[i] <- 1 - pchisq((n - 1) * (I - 1)*t.inertia * ((1 - sum(a[i, 
                3:M]^2))^(-1)) * (f[i, 1]^2/dmu[1, 1]^2 + f[i, 
                2]^2/dmu[2, 2]^2), df = (I - 1) * (J - 1))
        }
        else {
            pvalrow[i] <- 1 - pchisq((n - 1) * (I - 1)*t.inertia *  
                (f[i, 1]^2/dmu[1, 1]^2 + f[i, 2]^2/dmu[2, 2]^2), 
                df = (I - 1) * (J - 1))
        }
    }
    for (j in 1:J) {
        if (M > 2) {
       #     pvalcol[j] <- 1 - pchisq((n - 1) * (I - 1)*t.inertia * ((1/Jmass[j,j] - 
pvalcol[j] <- 1 - pchisq((n - 1) * (I - 1)*t.inertia * (1 - sum(b[j, 3:M]^2))^(-1) * (g[j, 1]^2/dmu[1, 1]^2 +g[j, 2]^2/dmu[2, 2]^2), df = (I - 1) * (J - 1))
        }
        else {
#            pvalcol[j] <- 1 - pchisq((n - 1) * (I - 1)*t.inertia * (1/Jmass[j,j]) * 
 pvalcol[j] <- 1 - pchisq((n - 1) * (I - 1)*t.inertia *(g[j, 1]^2/dmu[1, 1]^2 + g[j, 2]^2/dmu[2, 2]^2), df = (I - 1) * (J - 1))
#browser()
        }
    }
 #   summ.name <- c("HL Axis 1", "HL Axis 2", "Area", "P-value")
   summ.name <- c("HL Axis 1", "HL Axis 2", "Area")
#    row.summ <- cbind(hlax1.row, hlax2.row, area.row, pvalrow)
row.summ <- cbind(hlax1.row, hlax2.row, area.row) 
   dimnames(row.summ) <- list(paste(Inames[[1]]), paste(summ.name))
 #   col.summ <- cbind(hlax1.col, hlax2.col, area.col, pvalcol)
   col.summ <- cbind(hlax1.col, hlax2.col, area.col)

    dimnames(col.summ) <- list(paste(Jnames[[1]]), paste(summ.name))
  invisible  (list(eccentricity = eccentricity, row.summ = row.summ, 
        col.summ = col.summ))
 }
