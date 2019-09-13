plot.CAvariants<-function(x,  firstaxis=1, lastaxis=2, cex=0.8,cex.lab=0.8,prop=1,
plottype="biplot", biptype = "row",scaleplot=1,
posleg="topleft",pos=2,ell=FALSE,Mell=x$Mell,alpha=0.05, size=5,adj=c(0,0.5),...) {

## internal function to plot a  single picture
##

if ((firstaxis<1)|(firstaxis>x$r)) stop(paste("incorrect first axis =", firstaxis, "\n\n")) 
if (lastaxis>x$r) stop(paste("incorrect last axis =", lastaxis, "\n\n")) 
if (firstaxis>=lastaxis) stop(paste("last axis must be greater than first axis\n\n"))
#if (!any(plottype==c("classic","classical","c","biplot","bip","b"))) stop(paste("Must be specified the kind of graph: classic, or biplot"))

# Groups file must have no blank line at start and only one between sections
# group number   group name   symbol   colour   plot ellipse? 
n<-sum(x$Xtable)
I<-nrow(x$Xtable)
J<-ncol(x$Xtable)
rowgroup <- list(1:x$rows,rep(1,x$rows))
rowgrlab <- list(1,"","*","red","T")
colgroup <- list(1:x$cols,rep(1,x$cols))
colgrlab <- list(1,"","+","blue","T")
######################################################
# Plot row and col coordinates
#########################################
if ((plottype=="Classical")|(plottype=="classical")|(plottype=="classic")|(plottype=="c")) {
nthings<-x$cols
nvars<-x$rows
cord1<- x$Cprinccoord 
cord2<-x$Rprinccoord
#cord1<- x$Cstdcoord 
#cord2<-x$Rstdcoord
dmu=diag(x$inertias[,1])
inertiapc=round(x$inertias[,2],1) #inertia in percentage of row axes 
dimnames(cord1)[1]<-dimnames(x$Xtable)[2]
dimnames(cord2)[1]<-dimnames(x$Xtable)[1]
thinggroup<-colgroup
thinggrlab<-colgrlab
vargroup<-rowgroup
vargrlab<-rowgrlab
thinglabels<-x$collabels
varlabels<-x$rowlabels
main="Classical plot"
if ((x$catype=="DONSCA")|(x$catype=="DOCA")|(x$catype=="SOCA")|(x$catype=="SONSCA"))
{
cat("\n ERROR: NO CLASSICAL PLOT for ordered analysis. ONLY A BIPLOT can be constructed  (Please change 'plottype' and specify 'biptype')\n")
stop()
}
}#end classical plot
         if ((plottype=="Biplot")|(plottype=="biplot")|(plottype=="bip")|(plottype=="b")){
if ((biptype=="rows")|(biptype=="Rows")|(biptype=="row")|(biptype=="r")) 
{
plottype<-"biplot"
biptype<-"row"
cord1<-x$Rprinccoord*scaleplot
cord2<-x$Cstdcoord/scaleplot
nthings<-x$rows
nvars<-x$cols
thinglabels<-x$rowlabels
varlabels<-x$collabels
thinggroup<-rowgroup
thinggrlab<-rowgrlab
vargroup<-colgroup
vargrlab<-colgrlab
main<-"Row Isometric Biplot"
inertiapc=x$inertias[,2] #inertia of column poly
dmu=diag(x$inertias[,1])
dimnames(cord2)[1]<-dimnames(x$Xtable)[2]
dimnames(cord1)[1]<-dimnames(x$Xtable)[1]
if ((x$catype=="DONSCA")|(x$catype=="DOCA")|(x$catype=="SOCA")|(x$catype=="SONSCA"))
{
cord2<-x$Rprinccoord*scaleplot
cord1<-x$Cstdcoord/scaleplot
#cord2<- x$Rprinccoord
#cord1<-x$Cprinccoord
nthings<-x$cols
nvars<-x$rows
thinglabels<-x$collabels
varlabels<-x$rowlabels
thinggroup<-colgroup
thinggrlab<-colgrlab
vargroup<-rowgroup
vargrlab<-rowgrlab
inertiapc=round(x$inertias2[,2],1) #inertia of column poly
dmu=diag(x$inertias2[,1])
dimnames(cord2)[1]<-dimnames(x$Xtable)[1]
dimnames(cord1)[1]<-dimnames(x$Xtable)[2]
}#end catype

} #end bip row
if ((biptype=="cols")|(biptype=="Cols")|(biptype=="column")|(biptype=="col")) {
#plottype<-"biplot"
#biptype<-"columns"
#if (ell==TRUE) { scaleplot<-1}
if ((x$catype=="CA")|(x$catype=="NSCA")){
cord1<- x$Cprinccoord*scaleplot
cord2<-x$Rstdcoord/scaleplot
nthings<-x$cols
nvars<-x$rows
thinggroup<-colgroup
thinggrlab<-colgrlab
vargroup<-rowgroup
vargrlab<-rowgrlab
thinglabels<-x$collabels
varlabels<-x$rowlabels
main<-"Column Isometric Biplot"
inertiapc=round(x$inertias[,2],1) #inertia of row 
dmu=diag(x$inertias[,1])
dimnames(cord1)[1]<-dimnames(x$Xtable)[2]
dimnames(cord2)[1]<-dimnames(x$Xtable)[1]
}
if ((x$catype=="DONSCA")|(x$catype=="DOCA")|(x$catype=="SOCA")|(x$catype=="SONSCA"))
{
#if (ell==TRUE) {
#scaleplot<-1}
cord2<- x$Cprinccoord*scaleplot
cord1<-x$Rstdcoord/scaleplot
nthings<-x$rows
nvars<-x$cols
thinggroup<-rowgroup
thinggrlab<-rowgrlab
vargroup<-colgroup
vargrlab<-colgrlab
thinglabels<-x$rowlabels
varlabels<-x$collabels
inertiapc=round(x$inertias[,2],1) #inertia of row poly
dmu=diag(x$inertias[,1])
dimnames(cord1)[1]<-dimnames(x$Xtable)[1]
dimnames(cord2)[1]<-dimnames(x$Xtable)[2]
}#end catype
}#end bip column
}
###################################################################################ok without choice plottype
#repeat{
if ((x$catype=="DOCA")|(x$catype=="SOCA")|(x$catype=="SONSCA")|(x$catype=="DONSCA"))
{
 cat("\n Looking at the Trends of rows and columns\n")
#browser()
#trendplot(x@mj,(x@pcc), position="topleft",main="Original Trends of Centered #Column Profile",xlab="ordered scores")
#windows()
#trendplot(x@mi,t(x@pcc), position="topleft",main="Original Trends of Centered #Column Profile",xlab="ordered scores")
############################################################## reconstructed TREND
#windows()
plot.new()
trendplot(x$mj,(x$Trend), posleg=posleg,main="Reconstructed rows of the centred column profile",xlab="ordered scores",prop=prop)
dev.new()
trendplot(x$mi,t(x$Trend), posleg=posleg,main="Reconstructed columns of the centred column profile",xlab="ordered scores",prop=prop)
}
##############################################################
#library(scales)
#library (ggplot2)
#library(ggrepel)
#library(gridExtra)
categ<-NULL
frows <- data.frame(coord=cord1, labels=thinglabels, categ=rep("rows", nthings)) # build a dataframe to be used as input for plotting via ggplot2
  gcols <- data.frame(coord=cord2, labels=varlabels, categ=rep("cols", nvars)) # build a dataframe to be used as input for plotting via ggplot2
#-------------------------------------------------------------
  FGcord <- rbind(frows, gcols)                                       # build a dataframe to be used as input for plotting via   
   # vec <-cord2[, c(firstaxis, lastaxis)]
 #--------------------------------------------------
############################ 
plot.new()
if ((x$catype=="DONSCA")||(x$catype=="DOCA")||(x$catype=="SOCA")||(x$catype=="SONSCA"))
{
xmin <- min(FGcord[,firstaxis],FGcord[,lastaxis])
 xmax <- max(FGcord[,firstaxis],FGcord[,lastaxis])
ymin <- min(FGcord[,lastaxis],FGcord[,firstaxis])
ymax <- max(FGcord[,lastaxis],FGcord[,firstaxis])
 nv <- rep(0, nrow(cord1))
 CAplot <- ggplot(FGcord, aes(x=FGcord[,firstaxis], y=FGcord[,lastaxis]),type="b") + 
    geom_point(aes(color=categ, shape=categ), size=1.5) +
    geom_vline(xintercept = 0, linetype=2, color="gray") + 
    geom_hline(yintercept = 0, linetype=2, color="gray") + 
    labs(x=paste0("Principal Axis ",firstaxis,sep=" ", round(inertiapc[1],), "%"),y=paste0("Principal Axis ",lastaxis,sep= " ", round(inertiapc[2],0),"%"))  +  
    scale_x_continuous(limits = c(xmin, xmax)) +
    scale_y_continuous(limits = c(ymin, ymax)) + 
    theme(panel.background = element_rect(fill="white", colour="black")) + 
    scale_colour_manual(values=c("blue", "red")) + 
    coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE) + 
    geom_text_repel(data=FGcord, aes(colour=categ, label = labels), size = 3) +
geom_line(aes(group=categ,linetype=categ),lwd=.2)+
scale_linetype_manual(values=c("rows"="dashed","cols"="blank"))+
geom_segment(data=frows,aes(x=rep(0,c(nthings)),y=rep(0,c(nthings)),xend=frows[, firstaxis],yend=frows[, 
lastaxis]),colour=rep("blue",nthings))+
theme(legend.position="none")+
   ggtitle(" ") 
  grid.arrange(CAplot, ncol=1)
 }
######################
###############################
if ((plottype=="biplot")&&(x$catype=="CA")|(plottype=="biplot")&&(x$catype=="NSCA"))
{
xmin <- min(FGcord[,firstaxis],FGcord[,lastaxis])
 xmax <- max(FGcord[,firstaxis],FGcord[,lastaxis])
ymin <- min(FGcord[,lastaxis],FGcord[,firstaxis])
ymax <- max(FGcord[,lastaxis],FGcord[,firstaxis])
 nv <- rep(0, nrow(cord1))
 CAplot <- ggplot(FGcord, aes(x=FGcord[,firstaxis], y=FGcord[,lastaxis])) + 
    geom_point(aes(colour=categ, shape=categ), size=1.5) +
    geom_vline(xintercept = 0, linetype=2, color="gray") + 
    geom_hline(yintercept = 0, linetype=2, color="gray") + 
    labs(x=paste0("Principal Axis",firstaxis, sep=" ","", round(inertiapc[1],), "%"),  y=paste0("Principal Axis ",lastaxis,sep= " ", round(inertiapc[2],0),"%"))  +  
    scale_x_continuous(limits = c(xmin, xmax)) +
    scale_y_continuous(limits = c(ymin, ymax)) + 
    theme(panel.background = element_rect(fill="white", colour="black")) + 
    scale_color_manual(values=c("blue", "red")) + 
    coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE) + 
    geom_text_repel(data=FGcord, aes(colour=categ, label = labels), size = 3) +
geom_segment(data=frows,aes(x=rep(0,c(nthings)),y=rep(0,c(nthings)),xend=frows[, firstaxis],yend=frows[, 
lastaxis]),colour=rep("blue",nthings))+
theme(legend.position="none")+
    ggtitle(" ") 
  grid.arrange(CAplot, ncol=1)
}
##############################################################
if ((plottype=="classic")&&(x$catype=="CA")|(plottype=="classic")&&(x$catype=="NSCA"))
{
xmin <- min(FGcord[,firstaxis],FGcord[,lastaxis])
 xmax <- max(FGcord[,firstaxis],FGcord[,lastaxis])
ymin <- min(FGcord[,lastaxis],FGcord[,firstaxis])
ymax <- max(FGcord[,lastaxis],FGcord[,firstaxis])
 CAplot <- ggplot(FGcord, aes(x=FGcord[,firstaxis], y=FGcord[,lastaxis])) + 
    geom_point(aes(colour=categ, shape=categ), size=1.5) +
    geom_vline(xintercept = 0, linetype=2, color="gray") + 
    geom_hline(yintercept = 0, linetype=2, color="gray") + 
    labs(x=paste0("Principal Axis ", firstaxis,sep=" ", round(inertiapc[firstaxis]),"%"),  y=paste0("Principal Axis ", lastaxis,sep=" ", round(inertiapc[lastaxis]),"%" ))  +  
    scale_x_continuous(limits = c(xmin, xmax)) +
    scale_y_continuous(limits = c(ymin, ymax)) + 
    theme(panel.background = element_rect(fill="white", colour="black")) + 
    scale_color_manual(values=c("blue", "red")) + 
    coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE) + 
    geom_text_repel(data=FGcord, aes(colour=categ, label = labels), size = 3) +
theme(legend.position="none")+
    ggtitle(" ") 
  grid.arrange(CAplot, ncol=1)
}
#library("rgl")
#rgl::plot3d 
#if (plot3d==TRUE) {
#caplot3d(f=x$Rstdcoord,g=x$Cprinccoord,percIn=x$inertias[,2],size=size,adj=adj)
#}
##################
#cat("\nIncluding Beh's Confidence Ellipses\n")
################################################################################
if (ell==TRUE) {
if (((x$catype=="NSCA")|(x$catype=="CA")|(x$catype=="DOCA")|(x$catype=="SOCA")|(x$catype=="SONSCA")|(x$catype=="DONSCA")) & (plottype=="biplot")&(biptype=="row")|(biptype=="r")|(biptype=="rows")){
cordr<-cord2
cordc<-cord1
cord1<-cordr
cord2<-cordc
}
dev.new()
switch(x$catype, "CA"=caellipse(Xtable=x$Xtable,a1=firstaxis,a2=lastaxis,alpha=alpha,M=Mell,cex=cex,cex.lab=cex.lab,prop=prop,
Imass=x$Imass,Jmass=x$Jmass,a=x$Rstdcoord,b=x$Cstdcoord,g=cord1,fr=cord2,dmu=dmu,inertiapc=round(inertiapc,digits=1),
plottype=plottype,biptype=biptype,pos=pos,arrow=TRUE,length=0,graphy=TRUE,ell=TRUE), 

"SOCA"=caellipse(Xtable=x$Xtable,a1=firstaxis,a2=lastaxis,alpha=alpha,M=Mell,cex=cex,cex.lab=cex.lab,prop=prop,
Imass=x$Imass,Jmass=x$Jmass,a=solve(x$Imass^0.5)%*%x$Rstdcoord,b=solve(x$Jmass^0.5)%*%x$Cstdcoord,g=cord2,fr=cord1,dmu=dmu,inertiapc=round(inertiapc,digits=1),
plottype=plottype,biptype=biptype,pos=pos,arrow=FALSE,length=0,graphy=T,ell=TRUE),

"DOCA"=caellipse(Xtable=x$Xtable,a1=firstaxis,a2=lastaxis,alpha=alpha,M=Mell,cex=cex,cex.lab=cex.lab,prop=prop,
Imass=x$Imass,Jmass=x$Jmass,a=solve(x$Imass^0.5)%*%x$Rstdcoord,
b=solve(x$Jmass^0.5)%*%x$Cstdcoord,g=cord2,fr=cord1,dmu=dmu,inertiapc=round(inertiapc,digits=1),
plottype=plottype,biptype=biptype,pos=pos,arrow=FALSE,length=0,graphy=T,ell=TRUE), 

"NSCA"=nscaellipse(Xtable=x$Xtable,a1=firstaxis,a2=lastaxis,alpha=alpha,M=Mell,cex=cex,cex.lab=cex.lab,prop=prop,
Imass=x$Imass,Jmass=x$Jmass,a=x$Rstdcoord,
b=x$Cstdcoord,g=cord1,fr=cord2,dmu=dmu, tauden=x$tauden,
inertiapc=round(inertiapc,digits=1),
plottype=plottype,biptype=biptype,pos=pos,arrow=T,length=0,graphy=T,ell=TRUE),

"SONSCA"=nscaellipse(Xtable=x$Xtable,a1=firstaxis,a2=lastaxis,alpha=alpha,M= Mell,cex=cex,cex.lab=cex.lab,prop=prop,
Imass=x$Imass,Jmass=x$Jmass,a=x$Rstdcoord,
b=x$Cstdcoord,g=cord2,fr=cord1,dmu=dmu,tauden=x$tauden,inertiapc=round(inertiapc,digits=1),
plottype=plottype,biptype=biptype,pos=pos,arrow=FALSE,length=0,graphy=T,ell=TRUE), 

"DONSCA"=nscaellipse(Xtable=x$Xtable,a1=firstaxis,a2=lastaxis,alpha=alpha,M=Mell,cex=cex,cex.lab=cex.lab,prop=prop,
Imass=x$Imass,Jmass=x$Jmass,a=x$Rstdcoord,
b=x$Cstdcoord,g=cord2,fr=cord1,dmu=dmu,tauden=x$tauden,inertiapc=round(inertiapc,digits=1),
plottype=plottype,biptype=biptype,pos=pos,arrow=FALSE,length=0,graphy=T,ell=TRUE))

}#end if ellipse

}
