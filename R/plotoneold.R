plotone <-
function(a1,a2,plottype,things,nthings,nvars,Thingcoord,Varcoord,
                       inertiapc, thinggroup,thinggrlab,vargroup,vargrlab,thinglabels,varlabels,picsize,cex=0.8,cex.lab=0.8,
type="b" ,catype="CA" ,pos=2) {
#plot.new()
#plot(Thingcoord[ ,a1], Thingcoord[ ,a2], xlim=picsize, ylim=picsize, 
plot(0,0, xlim=picsize, ylim=picsize, 
   xlab=paste("Axis ", a1, "    ", inertiapc[a1], "%", sep=""), 
     ylab=paste("Axis ", a2, "    ", inertiapc[a2], "%", sep=""),  
        asp=1,pch="",cex=cex,cex.lab=cex.lab,type=type )
 #         asp=1, pch=thinggrlab[[3]][as.integer(thinggroup[[2]])], col=thinggrlab[[4]][as.integer(thinggroup[[2]])],cex=cex,cex.lab=cex.lab,type=type )
abline(h=0,v=0)
#if (as.integer(max(thinggroup[[2]]))>1) legend("topleft",thinggrlab[[2]],pch=thinggrlab[[3]])
if ((plottype=="Biplot")|(plottype=="biplot")|(plottype=="bip")) { 
  title(paste(things,"Biplot"))
  for (i in 1:nthings) { 
if ((catype=="NSCA")|(catype=="CA"))
{
lines(c(0,Thingcoord[i,a1]), c(0,Thingcoord[i,a2]), col=thinggrlab[[4]][as.integer(thinggroup[[2]][i])])}
#points(Thingcoord[i,a1],Thingcoord[i,a2],pch=thinggrlab[[3]][as.integer(thinggroup[[2]])],col=thinggrlab[[4]][as.integer(thinggroup[[2]])]) 
} 
 grat <- cbind(Thingcoord[,a1]/picsize[1],Thingcoord[,a1]/picsize[2],Thingcoord[,a2]/picsize[1],Thingcoord[,a2]/picsize[2],0.95)/0.95
 cl <- 1.05/apply(grat,1,max)
  text(cl*Thingcoord[ ,a1], cl*Thingcoord[ ,a2], labels=thinglabels, col=thinggrlab[[4]][as.integer(thinggroup[[2]])], pos=pos, cex=cex )
grat2 <- cbind(Varcoord[,a1]/picsize[1],Varcoord[,a1]/picsize[2],Varcoord[,a2]/picsize[1],Varcoord[,a2]/picsize[2],0.95)/0.95
  cl2 <- 1.05/apply(grat2,1,max)
points(Varcoord[ ,a1], Varcoord[ ,a2], asp=1,pch=vargrlab[[3]][as.integer(vargroup[[2]])], col=vargrlab[[4]][as.integer(vargroup[[2]])] ,cex=cex)
text(cl2*Varcoord[ ,a1], cl2*Varcoord[ ,a2], labels=varlabels, col=vargrlab[[4]][as.integer(vargroup[[2]])], pos=pos, cex=cex )
} else { # french
  title(paste(things,"Plot"))
 grat <- cbind(Thingcoord[,a1]/picsize[1],Thingcoord[,a1]/picsize[2],Thingcoord[,a2]/picsize[1],Thingcoord[,a2]/picsize[2],0.95)/0.95
 cl <- 1.05/apply(grat,1,max)
grat2 <- cbind(Varcoord[,a1]/picsize[1],Varcoord[,a1]/picsize[2],Varcoord[,a2]/picsize[1],Varcoord[,a2]/picsize[2],0.95)/0.95
  cl2 <- 1.05/apply(grat2,1,max)
  points(cl2*Varcoord[ ,a1],cl2* Varcoord[ ,a2], asp=1,  pch=vargrlab[[3]][as.integer(vargroup[[2]])], col=vargrlab[[4]][as.integer(vargroup[[2]])] ,cex=cex)
points(cl*Thingcoord[,a1],cl*Thingcoord[,a2],pch=thinggrlab[[3]][as.integer(thinggroup[[2]])],col=thinggrlab[[4]][as.integer(thinggroup[[2]])]) 
text(cl2*Varcoord[ ,a1],cl2* Varcoord[ ,a2], labels=varlabels, col=vargrlab[[4]][as.integer(vargroup[[2]])], pos=pos, cex=cex )
text(cl*Thingcoord[ ,a1], cl*Thingcoord[ ,a2], labels=thinglabels, col=thinggrlab[[4]][as.integer(thinggroup[[2]])], pos=pos, cex=cex )
#if (as.integer(max(vargroup[[2]]))>1) legend("topright",vargrlab[[2]],pch=vargrlab[[3]]) 
}

}
