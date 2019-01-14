nscaellipseprint <-function(Xtable, a1=1,a2=2,alpha = 0.05, M = 2,  
Imass,Jmass,a,b,f,g,dmu,tauden,inertiapc) 
{
  I <- nrow(Imass)
    J <- nrow(Jmass)
n<-sum(Xtable)
rowgroup <- list(1:I,rep(1,I))
rowgrlab <- list(1,"","*","blue","T")
colgroup <- list(1:J,rep(1,J))
colgrlab <- list(1,"","+","red","T")
    Inames <- dimnames(f)[1]
    Jnames <- dimnames(g)[1]
     t.inertia <- (sum(dmu)/tauden)*(n-1)*(I-1)
   #    t.inertia2 <- (sum(dmu^2))
 dmu<-sqrt(dmu)
#   chisq.val <- (qchisq(1 - alpha, df = (I - 1) * (J - 1)))*tauden/((n - 1) * (I - 1))
   chisq.val <- (qchisq(1 - alpha, df = (I - 1) * (J - 1))) 
   hlax1.row <- vector(mode = "numeric", length = I)
    hlax2.row <- vector(mode = "numeric", length = I)
    hlax1.col <- vector(mode = "numeric", length = J)
    hlax2.col <- vector(mode = "numeric", length = J)
if (M > 2) {
        for (i in 1:I) {
            hlax1.row[i] <- dmu[a1, a1] * sqrt(abs((chisq.val/t.inertia) * 
                (1 - sum(a[i, 3:M]^2))))
            hlax2.row[i] <- dmu[a2, a2] * sqrt(abs((chisq.val/t.inertia) * 
                (1 - sum(a[i, 3:M]^2))))
        }
        for (j in 1:J) {
            hlax1.col[j] <- dmu[a1, a1] * sqrt(abs((chisq.val/t.inertia) * 
                (1/Jmass[j,j] - sum(b[j, 3:M]^2))))
            hlax2.col[j] <- dmu[a2, a2] * sqrt(abs((chisq.val/t.inertia) * 
                (1/Jmass[j,j] - sum(b[j, 3:M]^2))))
        }
    }
    else {
        for (i in 1:I) {
            hlax1.row[i] <- dmu[a1, a1] * sqrt(chisq.val/t.inertia)
            hlax2.row[i] <- dmu[a2, a2] * sqrt(chisq.val/t.inertia) 
        }
        for (j in 1:J) {
            hlax1.col[j] <- dmu[a1, a1] * sqrt((chisq.val/t.inertia) * 
                1/Jmass[j,j])
            hlax2.col[j] <- dmu[a2, a2] * sqrt((chisq.val/t.inertia) * 
                1/Jmass[j,j])
        }
    }

#----------------------------------------
   eccentricity <- abs(1 - (dmu[a2, a2]/dmu[a1, a1])^2)^(1/2)
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
            pvalrow[i] <- 1 - pchisq(t.inertia * ((1 - sum(a[i, 
                3:M]^2))^(-1)) * (f[i, a1]^2/dmu[a1, a1]^2 + f[i, 
                a2]^2/dmu[a2, a2]^2), df =  (J - 1))
        }
        else {
            pvalrow[i] <- 1 - pchisq(t.inertia *
                (f[i, a1]^2/dmu[a1, a1]^2 + f[i, a2]^2/dmu[a2, a2]^2), 
                df =  (J - 1))
        }
    }
    for (j in 1:J) {
        if (M > 2) {
    pvalcol[j] <- 1 - pchisq(t.inertia * (1/Jmass[j,j] - sum(b[j, 3:M]^2))^(-1) * (g[j, a1]^2/dmu[a1, a1]^2 +g[j, a2]^2/dmu[a2,a2]^2), df = (I - 1) )
        }
        else {
  pvalcol[j] <- 1 - pchisq(t.inertia* (Jmass[j,j])*(g[j, a1]^2/dmu[a1, a1]^2 + g[j, a2]^2/dmu[a2, a2]^2), df = (I - 1) )
#browser()
        }
    }
    summ.name <- c("HL Axis 1", "HL Axis 2", "Area", "P-value")
    row.summ <- cbind(hlax1.row, hlax2.row, area.row, pvalrow)
   dimnames(row.summ) <- list(paste(Inames[[1]]), paste(summ.name))
    col.summ <- cbind(hlax1.col, hlax2.col, area.col, pvalcol)
    dimnames(col.summ) <- list(paste(Jnames[[1]]), paste(summ.name))
invisible(  list(eccentricity = eccentricity, row.summ = row.summ, 
        col.summ = col.summ))
 }
