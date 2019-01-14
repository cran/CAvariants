trendplot <-
function(mf,mg, cex = 1,  cex.lab = 0.8,prop=0.5,posleg="topleft",
xlab="First Principal Axis",ylab="Second Principal Axis")
{
#------------------------------------------------------------------------------
# f and g are the row and column coordinates
#------------------------------------------------------------------------------
nrows<-dim(mg)[[1]]
ncols<-dim(mg)[[2]]
leg.txt<-dimnames(mg)[[1]]
colsymb<-c(1:nrows)
gt<-t(mg)
plot(mf,mg[1,],type="b", ylim = range(gt[1:ncols,],mg[1:nrows,])/prop, xlab = xlab, ylab = ylab, cex = cex, cex.lab = cex.lab, col=1)
abline(h=0,lty=3)
   for (i in 1:(nrows)){
lines(mf,mg[i,],type="b",pch=i,  col=i)
}
legend(x=posleg,legend=leg.txt,col=colsymb,pch=c(1:nrows),bty="o",cex=.8)

}
