trendplot <-
function(f,g, cex = 1,  cex.lab = 0.8, main=" ",prop=0.5,posleg="right",xlab="First Axis",ylab="Second Axis")
{
#------------------------------------------------------------------------------
# f and g are the row and column coordinates
#------------------------------------------------------------------------------
par(mar=c(5,4,4,8),xpd=TRUE)
nrows<-dim(g)[[1]]
ncols<-dim(g)[[2]]
leg.txt<-dimnames(g)[[2]]
leg.txt1 <-dimnames(g)[[1]]
colsymb<-c(1:nrows)
gt<-t(g)
plot(f,g[1,],type="b", ylim = range(gt[1:ncols,],g[1:nrows,])/0.5, xlab = xlab, ylab = ylab, cex = cex, cex.lab = cex.lab, main=main, col=1,xaxt="n")
axis(1, at= 1:ncols, labels=leg.txt)
#axis(1, at labels=leg.txt)
#abline(h=0,lty=3)
   for (i in 1:(nrows)){
lines(f,g[i,],type="b",pch=i,  col=i)
}
legend(x=posleg,legend=leg.txt1, inset=c(-0.4,0),col=c(1:nrows),pch=c(1:nrows),bty="o",cex=.8)
#----------------------------------------------------------------------------------------
   }



