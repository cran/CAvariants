plot.CAvariants<-function(x,  firstaxis=1, lastaxis=2, thirdaxis=3, cex=0.8,cex.lab=0.8,
plottype="biplot", biptype = "row",scaleplot=1,
posleg="topleft",pos=2,ell=FALSE,alpha=0.05,plot3d =FALSE,size1=1.5,size2=3,invproj=TRUE,...) {
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
######################################################
# Plot row and col coordinates
#########################################
if ((plottype=="Classical")|(plottype=="classical")|(plottype=="classic")|(plottype=="Classic")) {
nthings<-x$cols
nvars<-x$rows
cord1<- x$Cprinccoord*scaleplot 
cord2<-x$Rprinccoord/scaleplot
#cord1<- x$Cstdcoord 
#cord2<-x$Rstdcoord
dmu=diag(x$inertias[,1])
inertiapc=round(x$inertias[,2],1) #inertia in percentage of row axes 
dimnames(cord1)[1]<-dimnames(x$Xtable)[2]
dimnames(cord2)[1]<-dimnames(x$Xtable)[1]
thinglabels<-x$collabels
varlabels<-x$rowlabels
main="Classical plot"
if ((x$catype=="DONSCA")|(x$catype=="DOCA")|(x$catype=="SOCA")|(x$catype=="SONSCA"))
{
cat("\n ERROR: NO CLASSICAL PLOT for ordered analysis. ONLY A BIPLOT can be constructed  (Please change 'plottype' and specify 'biptype')\n")
stop()
}
#---------------------------------------------------------------------------------------
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
main<-"Row Isometric Biplot"
inertiapc=x$inertias[,2] #inertia of rows
dmu=diag(x$inertias[,1])
dimnames(cord2)[1]<-dimnames(x$Xtable)[2]
dimnames(cord1)[1]<-dimnames(x$Xtable)[1]
if ((x$catype=="DONSCA")|(x$catype=="DOCA")|(x$catype=="SOCA")|(x$catype=="SONSCA"))
{
cord2<-x$Rprinccoord*scaleplot
cord1<-x$Cstdcoord/scaleplot
nthings<-x$cols
nvars<-x$rows
thinglabels<-x$collabels
varlabels<-x$rowlabels
inertiapc=round(x$inertias2[,2],1) #inertia of rows
dmu=diag(x$inertias2[,1])
dimnames(cord2)[1]<-dimnames(x$Xtable)[1]
dimnames(cord1)[1]<-dimnames(x$Xtable)[2]
}#end catype

} #end bip row
if ((biptype=="cols")|(biptype=="Cols")|(biptype=="column")|(biptype=="col")) {
if ((x$catype=="CA")|(x$catype=="NSCA")){
cord1<- x$Cprinccoord*scaleplot
cord2<-x$Rstdcoord/scaleplot
nthings<-x$cols
nvars<-x$rows
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
cord2<- x$Cprinccoord*scaleplot
cord1<-x$Rstdcoord/scaleplot
nthings<-x$rows
nvars<-x$cols
thinglabels<-x$rowlabels
varlabels<-x$collabels
inertiapc=round(x$inertias[,2],1) #inertia of cols
dmu=diag(x$inertias[,1])
dimnames(cord1)[1]<-dimnames(x$Xtable)[1]
dimnames(cord2)[1]<-dimnames(x$Xtable)[2]
}#end catype
}#end bip column
}
###################################################################################ok without choice plottype
#if ((x$catype=="DOCA")|(x$catype=="SOCA")|(x$catype=="SONSCA")|(x$catype=="DONSCA"))
#{
# cat("\n Looking at the Trends of rows and columns\n")
############################################################## reconstructed TREND
#trendplot(x$mj,(x$Trend), posleg=posleg, xlab="ordered scores",prop=prop)
#dev.new()
#trendplot(x$mi,t(x$Trend), posleg=posleg,xlab="ordered scores",prop=prop)
#}
##############################################################
##################
#library(scales)
#library (ggplot2)
#library(ggrepel)
#library(gridExtra)
#---------------------------------------------------------------------------------
categ<-NULL
frows <- data.frame(coord=cord1, labels=thinglabels, categ=rep("rows", nthings)) # build a dataframe to be used as input for plotting via ggplot2
gcols <- data.frame(coord=cord2, labels=varlabels, categ=rep("cols", nvars)) # build a dataframe to be used as input for plotting via ggplot2
FGcord <- rbind(frows, gcols)                                       # build a dataframe to be used as input for plotting via   
############################################################ 
if (((x$catype=="DONSCA")||(x$catype=="DOCA"))&&((biptype=="column")&(plottype=="biplot")))
{
caplotord(frows=frows,gcols=gcols,firstaxis=firstaxis,lastaxis=lastaxis,nseg=nvars,inertiapc=inertiapc,thingseg=gcols,col1="red",
col2="blue",col3="blue",size1=size1,size2=size2)
if (invproj==TRUE){
caplotord(frows=frows,gcols=gcols,firstaxis=firstaxis,lastaxis=lastaxis,nseg=nthings,inertiapc=inertiapc,thingseg=frows,col1="red",
col2="blue",col3="red",size1=size1,size2=size2)
}
 }#end catype
if (((x$catype=="SONSCA")||(x$catype=="SOCA"))&&((biptype=="column")&(plottype=="biplot")))
{
caplot(frows=frows,gcols=gcols,firstaxis=firstaxis,lastaxis=lastaxis,nseg=nvars,inertiapc=inertiapc,thingseg=gcols,col1="red",col2="blue",
col3="blue",size1=size1,size2=size2)
if (invproj==FALSE){
caplot(frows=frows,gcols=gcols,firstaxis=firstaxis,lastaxis=lastaxis,nseg=nthings,inertiapc=inertiapc,thingseg=frows,col1="red",col2="blue",
col3="red",size1=size1,size2=size2)
}
 }
###############################################################
if (((x$catype=="DONSCA")||(x$catype=="DOCA")||(x$catype=="SOCA")||(x$catype=="SONSCA"))&&((biptype=="row")&(plottype=="biplot")))
{
caplotord(frows=gcols,gcols=frows,firstaxis=firstaxis,lastaxis=lastaxis,nseg=nvars,inertiapc=inertiapc,thingseg=gcols,col1="red",
col2="blue",col3="red",size1=size1,size2=size2)
 
if (invproj==FALSE){
caplotord(frows=gcols,gcols=frows,firstaxis=firstaxis,lastaxis=lastaxis,nseg=nthings,inertiapc=inertiapc,thingseg=frows,col1="red",
col2="blue",col3="blue",size1=size1,size2=size2)
}
}
#-----------------------------------------------------------
if (((x$catype=="NSCA")||(x$catype=="CA"))&&((biptype=="row")&(plottype=="biplot"))) 
{
caplot(frows=gcols,gcols=frows,firstaxis=firstaxis,lastaxis=lastaxis,nseg=nthings,inertiapc=inertiapc,thingseg=frows,col1="blue",col2="red",
col3="red",size1=size1,size2=size2)
if (invproj==FALSE){
caplot(frows=gcols,gcols=frows,firstaxis=firstaxis,lastaxis=lastaxis,nseg=nvars,inertiapc=inertiapc,thingseg=gcols,col1="blue",col2="red",
col3="blue",size1=size1,size2=size2)
}
  }
###############################
if (((x$catype=="NSCA")||(x$catype=="CA"))&&((biptype=="column")&(plottype=="biplot"))) 
{
caplot(frows=frows,gcols=gcols,firstaxis=firstaxis,lastaxis=lastaxis,nseg=nthings,inertiapc=inertiapc,thingseg=frows,col1="blue",col2="red",
col3="blue",size1=size1,size2=size2)
if (invproj==FALSE){
caplot(frows=frows,gcols=gcols,firstaxis=firstaxis,lastaxis=lastaxis,nseg=nvars,inertiapc=inertiapc,thingseg=gcols,col1="blue",col2="red",
col3="red",size1=size1,size2=size2)
}
}
##############################################################
if ((plottype=="classic")&&(x$catype=="CA")|(plottype=="classic")&&(x$catype=="NSCA"))
{
categ<-NULL
frows <- data.frame(coord=cord1, labels=thinglabels, categ=rep("rows", nthings)) # build a dataframe to be used as input for plotting via ggplot2
gcols <- data.frame(coord=cord2, labels=varlabels, categ=rep("cols", nvars)) # build a dataframe to be used as input for plotting via ggplot2
#-------------------------------------------------------------
FGcord <- rbind(frows, gcols)                                       # build a dataframe to be used as input for plotting via   
xmin <- min(FGcord[,firstaxis],FGcord[,lastaxis])
 xmax <- max(FGcord[,firstaxis],FGcord[,lastaxis])
ymin <- min(FGcord[,lastaxis],FGcord[,firstaxis])
ymax <- max(FGcord[,lastaxis],FGcord[,firstaxis])

 CAplot <- ggplot(FGcord, aes(x=FGcord[,firstaxis], y=FGcord[,lastaxis])) + 
    geom_point(aes(colour=categ, shape=categ), size=size1) +
    geom_vline(xintercept = 0, linetype=2, color="gray") + 
    geom_hline(yintercept = 0, linetype=2, color="gray") + 
    labs(x=paste0("Dimension ", firstaxis,sep=" (", round(inertiapc[firstaxis],1),"%) "),  y=paste0("Dimension ", lastaxis,sep=" (", round(inertiapc[lastaxis],1),"%)" ))  +  
    scale_x_continuous(limits = c(xmin, xmax)) +
    scale_y_continuous(limits = c(ymin, ymax)) + 
    theme(panel.background = element_rect(fill="white", colour="black")) + 
    scale_color_manual(values=c("blue", "red")) + 
    coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE) + 
    geom_text_repel(data=FGcord, aes(colour=categ, label = labels), size = size2) +
theme(legend.position="none")+
    ggtitle(" ") 
  grid.arrange(CAplot, ncol=1)
}
#cat("\nIncluding Beh's Confidence Ellipses\n")
################################################################################
if (ell==TRUE) {
cord1<-x$Cprinccoord*scaleplot #check here!!
cord2<-x$Rprinccoord/scaleplot
#if (((x$catype=="DOCA")|(x$catype=="SOCA")|(x$catype=="SONSCA")|(x$catype=="DONSCA")) & (plottype=="biplot")&(biptype=="row")|(biptype=="r")|(biptype=="rows")){
if ((x$catype=="DOCA")|(x$catype=="SOCA")|(x$catype=="SONSCA")|(x$catype=="DONSCA")){
cordr<-cord2
cordc<-cord1
cord1<-cordr
cord2<-cordc
}
#dev.new()
vcaellipse(t.inertia=x$t.inertia,inertias=x$inertias[,1],inertiapc=x$inertias[,2],cord1=x$Rprinccoord,cord2=x$Cprinccoord,a=x$Rstdcoord,b=x$Cstdcoord,firstaxis=firstaxis,lastaxis=lastaxis,n=x$n,M=x$M,Imass=x$Imass,Jmass=x$Jmass) 

}#end if ellipse

#library("plot3D")
if (plot3d==TRUE) {
#coordR<-cord1
#coordC<-cord2
inertiaper=x$inertias[,2]
 caplot3d(coordR=cord1,coordC =cord2,inertiaper=x$inertias[,2])
}

}