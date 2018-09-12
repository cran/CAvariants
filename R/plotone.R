plotone <-
function(a1,a2,plottype,things,nthings,nvars,Thingcoord,Varcoord,
                       inertiapc,
                       thinggroup,thinggrlab,vargroup,vargrlab,thinglabels,varlabels,picsize,cex=0.8,cex.lab=0.8,type ,catype ,pos=2) {
 grat <- cbind(Thingcoord[,a1]/picsize[1],Thingcoord[,a1]/picsize[2],Thingcoord[,a2]/picsize[1],Thingcoord[,a2]/picsize[2],0.95)/0.95
 cl <- 1.05/apply(grat,1,max)
grat2 <- cbind(Varcoord[,a1]/picsize[1],Varcoord[,a1]/picsize[2],Varcoord[,a2]/picsize[1],Varcoord[,a2]/picsize[2],0.95)/0.95
  cl2 <- 1.05/apply(grat2,1,max)
#browser()
plot(cl*Thingcoord[ ,a1], cl*Thingcoord[ ,a2], xlim=picsize, ylim=picsize, 
     xlab=paste("Axis ", a1, "    ", inertiapc[a1], "%", sep=""), 
     ylab=paste("Axis ", a2, "    ", inertiapc[a2], "%", sep=""),  
          asp=1, pch=thinggrlab[[3]][as.integer(thinggroup[[2]])], col=thinggrlab[[4]][as.integer(thinggroup[[2]])],cex=cex,cex.lab=cex.lab,type=type,lty=3 )
abline(h=0,v=0)
if ((plottype=="Biplot")|(plottype=="biplot")|(plottype=="bip")|(plottype=="b")) { 
  title(paste(things,"Biplot"))
 text(cl*Thingcoord[ ,a1], cl*Thingcoord[ ,a2], labels=thinglabels, col=thinggrlab[[4]][as.integer(thinggroup[[2]])], pos=pos, cex=cex )
for (i in 1:nvars) {
#points(cl2*Varcoord[ ,a1], cl2*Varcoord[ ,a2], asp=1,pch=vargrlab[[3]][as.integer(vargroup[[2]])], col=vargrlab[[4]][as.integer(vargroup[[2]])] ,cex=cex)
lines(c(0,cl2*Varcoord[i,a1]), c(0,cl2*Varcoord[i,a2]), col=vargrlab[[4]][as.integer(vargroup[[2]][i])])
}
text(cl2*Varcoord[ ,a1], cl2*Varcoord[ ,a2], labels=varlabels, col=vargrlab[[4]][as.integer(vargroup[[2]])], pos=pos, cex=cex )
} else { # french
  title(paste(things,"Plot"))
points(cl2*Varcoord[ ,a1],cl2* Varcoord[ ,a2], asp=1, 
         pch=vargrlab[[3]][as.integer(vargroup[[2]])], col=vargrlab[[4]][as.integer(vargroup[[2]])] ,cex=cex)
text(cl2*Varcoord[ ,a1], cl2*Varcoord[ ,a2], labels=varlabels, col=vargrlab[[4]][as.integer(vargroup[[2]])], pos=pos, cex=cex )
text(cl*Thingcoord[ ,a1],cl* Thingcoord[ ,a2], labels=thinglabels, col=thinggrlab[[4]][as.integer(thinggroup[[2]])], pos=pos, cex=cex )
}
}
