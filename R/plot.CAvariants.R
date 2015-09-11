plot.CAvariants<-function(x,  firstaxis=1, lastaxis=2, cex=1,prop=1,
plottype="biplot", biptype = "rows",M=2,scaleplot=3,
posleg="topleft",pos=2,ell=FALSE,...) {

## internal function to plot a  single picture
##

if ((firstaxis<1)|(firstaxis>x$cacorpo@maxaxes-1)) stop(paste("incorrect first axis =", firstaxis, "\n\n")) 
if (lastaxis>x$cacorpo@maxaxes) stop(paste("incorrect last axis =", lastaxis, "\n\n")) 
if (firstaxis>=lastaxis) stop(paste("last axis must be greater than first axis\n\n"))
#if (!any(plottype==c("classic","classical","c","biplot","bip","b"))) stop(paste("Must be specified the kind of graph: classic, or biplot"))

# Groups file must have no blank line at start and only one between sections
# group number   group name   symbol   colour   plot ellipse? 

rowgroup <- list(1:x$cacorpo@rows,rep(1,x$cacorpo@rows))
rowgrlab <- list(1,"","*","red","T")
colgroup <- list(1:x$cacorpo@cols,rep(1,x$cacorpo@cols))
colgrlab <- list(1,"","+","blue","T")
######################################################
# Plot row and col coordinates
#########################################
if ((plottype=="Classical")|(plottype=="classical")|(plottype=="classic")|(plottype=="c")) {
nthings<-x$cacorpo@cols
nvars<-x$cacorpo@rows
cord1<- x$cacorpo@Cprinccoord 
cord2<-x$cacorpo@Rprinccoord
inertiapc=x$cacorpo@inertias[,2] #inertia in percentage of row axes 
dimnames(cord1)[1]<-dimnames(x$cacorpo@DataMatrix)[2]
dimnames(cord2)[1]<-dimnames(x$cacorpo@DataMatrix)[1]
thinggroup<-colgroup
thinggrlab<-colgrlab
vargroup<-rowgroup
vargrlab<-rowgrlab
thinglabels<-x$cacorpo@collabels
varlabels<-x$cacorpo@rowlabels
main="Classical plot"
if ((x$cacorpo@catype=="DONSCA")|(x$cacorpo@catype=="DOCA")|(x$cacorpo@catype=="SOCA")|(x$cacorpo@catype=="SONSCA"))
{
cat("\n NO CLASSICAL PLOT ONLY BIPLOT (coordinates should be biplot coordinates, please change plottype)\n")
cord1<-x$cacorpo@Cstdcoord/scaleplot
cord2<-x$cacorpo@Rprinccoord*scaleplot
inertiapc=x$cacorpo@inertias2[,2] #inertia in percentage of row axes 
}

}

         if ((plottype=="Biplot")|(plottype=="biplot")|(plottype=="bip")|(plottype=="b")){
if ((biptype=="rows")|(biptype=="Rows")|(biptype=="row")|(biptype=="r")) 
{
plottype<-"biplot"
#biptype<-"row"
cord1<-x$cacorpo@Rprinccoord*scaleplot
cord2<-x$cacorpo@Cstdcoord/scaleplot
nthings<-x$cacorpo@rows
nvars<-x$cacorpo@cols
thinglabels<-x$cacorpo@rowlabels
varlabels<-x$cacorpo@collabels
thinggroup<-rowgroup
thinggrlab<-rowgrlab
vargroup<-colgroup
vargrlab<-colgrlab
main<-"Row Isometric Biplot"
inertiapc=x$cacorpo@inertias[,2] #inertia of column poly
dimnames(cord2)[1]<-dimnames(x$cacorpo@DataMatrix)[2]
dimnames(cord1)[1]<-dimnames(x$cacorpo@DataMatrix)[1]

if ((x$cacorpo@catype=="DONSCA")|(x$cacorpo@catype=="DOCA")|(x$cacorpo@catype=="SOCA")|(x$cacorpo@catype=="SONSCA"))
{
cord2<-x$cacorpo@Rprinccoord*scaleplot
cord1<-x$cacorpo@Cstdcoord/scaleplot
nthings<-x$cacorpo@cols
nvars<-x$cacorpo@rows
thinglabels<-x$cacorpo@collabels
varlabels<-x$cacorpo@rowlabels
thinggroup<-colgroup
thinggrlab<-colgrlab
vargroup<-rowgroup
vargrlab<-rowgrlab
inertiapc=x$cacorpo@inertias2[,2] #inertia of column poly
dimnames(cord2)[1]<-dimnames(x$cacorpo@DataMatrix)[1]
dimnames(cord1)[1]<-dimnames(x$cacorpo@DataMatrix)[2]

}#end catype

} #end bip row
else{
plottype<-"biplot"
#biptype<-"columns"
cord1<- x$cacorpo@Cprinccoord*scaleplot
cord2<-x$cacorpo@Rstdcoord/scaleplot
nthings<-x$cacorpo@cols
nvars<-x$cacorpo@rows
thinggroup<-colgroup
thinggrlab<-colgrlab
vargroup<-colgroup
vargrlab<-rowgrlab
thinglabels<-x$cacorpo@collabels
varlabels<-x$cacorpo@rowlabels
main<-"Column Isometric Biplot"
inertiapc=x$cacorpo@inertias[,2] #inertia of row 
dimnames(cord1)[1]<-dimnames(x$cacorpo@DataMatrix)[2]
dimnames(cord2)[1]<-dimnames(x$cacorpo@DataMatrix)[1]

if ((x$cacorpo@catype=="DONSCA")|(x$cacorpo@catype=="DOCA")|(x$cacorpo@catype=="SOCA")|(x$cacorpo@catype=="SONSCA"))
{
cord2<- x$cacorpo@Cprinccoord*scaleplot
cord1<-x$cacorpo@Rstdcoord/scaleplot
nthings<-x$cacorpo@rows
nvars<-x$cacorpo@cols
thinggroup<-rowgroup
thinggrlab<-rowgrlab
vargroup<-colgroup
vargrlab<-colgrlab
thinglabels<-x$cacorpo@rowlabels
varlabels<-x$cacorpo@collabels

inertiapc=x$cacorpo@inertias[,2] #inertia of row poly
dimnames(cord1)[1]<-dimnames(x$cacorpo@DataMatrix)[1]
dimnames(cord2)[1]<-dimnames(x$cacorpo@DataMatrix)[2]

}#end catype

}#end bip column
}


###################################################################################ok without choice plottype
#repeat{

if ((x$cacorpo@catype=="DOCA")|(x$cacorpo@catype=="SOCA")|(x$cacorpo@catype=="SONSCA")|(x$cacorpo@catype=="DONSCA"))
{
 cat("\n Looking at the Trends of rows and columns\n")
#browser()
#trendplot(x@mj,(x@pcc), position="topleft",main="Original Trends of Centered #Column Profile",xlab="ordered scores")
#windows()
#trendplot(x@mi,t(x@pcc), position="topleft",main="Original Trends of Centered #Column Profile",xlab="ordered scores")
############################################################## reconstructed TREND
#windows()
trendplot(x$cacorpo@mj,(x$cacorpo@Trend), posleg=posleg,main="Reconstructed rows of the centred column profile",xlab="ordered scores",prop=prop)
dev.new()
trendplot(x$cacorpo@mi,t(x$cacorpo@Trend), posleg=posleg,main="Reconstructed columns of the centred column profile",xlab="ordered scores",prop=prop)
#browser()
}

##############################################
picsize1<-c(range(cord1[,c(firstaxis,lastaxis)], cord2[,c(firstaxis,lastaxis)])/prop)
picsizeR<-c(range(cord1[,c(firstaxis,lastaxis)], cord1[,c(firstaxis,lastaxis)])/prop)
picsize2<-c(range(cord2[,c(firstaxis,lastaxis)], cord2[,c(firstaxis,lastaxis)])/prop)

if (picsize1[1]>=picsize1[2]) stop(paste("incorrect axis scale picsize =", picsize1[1], picsize1[2], "\n\n"))
########################################################################################## 
if ((x$cacorpo@catype=="DONSCA")||(x$cacorpo@catype=="DOCA"))
{
plotone (firstaxis,lastaxis,plottype=plottype,things=x$cacorpo@catype,nthings,nvars,cord1,cord2,
inertiapc=inertiapc, thinggroup,thinggrlab,vargroup,vargrlab, thinglabels, varlabels,picsize=picsize1,cex=cex,type="b",
catype=x$cacorpo@catype,pos=pos) 
}
######################
if ((x$cacorpo@catype=="SOCA")||(x$cacorpo@catype=="SONSCA"))
{
if (biptype=="row"){type="b"}
else {type="p"}
plotone (firstaxis,lastaxis,plottype=plottype,things=x$cacorpo@catype,nthings,nvars,cord1,cord2,
inertiapc=inertiapc, thinggroup,thinggrlab,vargroup,vargrlab, thinglabels, varlabels,picsize=picsize1,cex=cex,catype=x$cacorpo@catype,type=type,pos=pos) 
}

###############################
if ((x$cacorpo@catype=="CA")||(x$cacorpo@catype=="NSCA"))
{
plotone (firstaxis,lastaxis,plottype=plottype,things=x$cacorpo@catype,nthings,nvars,cord1,cord2,
inertiapc=inertiapc, thinggroup,thinggrlab,vargroup,vargrlab, thinglabels, varlabels,picsize=picsize1,cex=cex,type="p",catype=x$cacorpo@catype,pos=pos) 
}
################################################################################
#cat("\nIncluding Beh's Confidence Ellipses\n")
#####################################################
if (ell==TRUE) {
if ((x$cacorpo@catype=="CA")&(plottype=="biplot")&(biptype=="row")){
cordr<-cord2
cordc<-cord1
cord1<-cordr
cord2<-cordc
}
if ((x$cacorpo@catype=="NSCA")&(plottype=="biplot")&(biptype=="row")){
cordr<-cord2
cordc<-cord1
cord1<-cordr
cord2<-cordc
}

if (((x$cacorpo@catype=="DOCA")|(x$cacorpo@catype=="SOCA")|(x$cacorpo@catype=="SONSCA")|(x$cacorpo@catype=="DONSCA")) & (plottype=="biplot")&((biptype=="column")|(biptype=="col"))){
cordr<-cord2
cordc<-cord1
cord1<-cordr
cord2<-cordc
}
plot.new()

switch(x$cacorpo@catype, "CA"=caellipse(N=x$cacorpo@DataMatrix,a1=firstaxis,a2=lastaxis,M=M,prop=prop,
Imass=x$cacorpo@Imass,Jmass=x$cacorpo@Jmass,a=x$cacorpo@S@Caxes,b=x$cacorpo@S@Raxes,g=cord1,f=cord2,dmu=diag(x$cacorpo@inertias[,1]),inertiapc=x$cacorpo@inertias[,2],
plottype=plottype,biptype=biptype,pos=pos,arrow=T), 

"SOCA"=caellipse(N=x$cacorpo@DataMatrix,a1=firstaxis,a2=lastaxis,M=M,prop=prop,
Imass=x$cacorpo@Imass,Jmass=x$cacorpo@Jmass,a=x$cacorpo@S@Caxes,b=x$cacorpo@S@Raxes,g=cord2,f=cord1,dmu=diag(x$cacorpo@inertias2[,1]),inertiapc=x$cacorpo@inertias2[,2],
plottype=plottype,biptype=biptype,pos=pos,arrow=F), 

"DOCA"=caellipse(N=x$cacorpo@DataMatrix,a1=firstaxis,a2=lastaxis,M=M,prop=prop,
Imass=x$cacorpo@Imass,Jmass=x$cacorpo@Jmass,a=x$cacorpo@S@Caxes,
b=x$cacorpo@S@Raxes,g=cord2,f=cord1,dmu=diag(x$cacorpo@inertias2[,1]),inertiapc=x$cacorpo@inertias2[,2],
plottype=plottype,biptype=biptype,pos=pos,arrow=F), 

"NSCA"=nsca.ellipse(x$cacorpo@DataMatrix,a1=firstaxis,a2=lastaxis,M=M,prop=prop,
Imass=x$cacorpo@Imass,Jmass=x$cacorpo@Jmass,a=x$cacorpo@S@Caxes,
b=x$cacorpo@S@Raxes,g=cord1,f=cord2,dmu=diag(x$cacorpo@inertias[,1]),
inertiapc=x$cacorpo@inertias[,2],
plottype=plottype,biptype=biptype,pos=pos,arrow=T),

"SONSCA"=nsca.ellipse(x$cacorpo@DataMatrix,a1=firstaxis,a2=lastaxis,M= M,prop=prop,
Imass=x$cacorpo@Imass,Jmass=x$cacorpo@Jmass,a=x$cacorpo@S@Caxes,
b=x$cacorpo@S@Raxes,g=cord2,f=cord1,dmu=diag(x$cacorpo@inertias2[,1]),inertiapc=x$cacorpo@inertias2[,2],
plottype=plottype,biptype=biptype,pos=pos,arrow=F), 

"DONSCA"=nsca.ellipse(x$cacorpo@DataMatrix,a1=firstaxis,a2=lastaxis,M=M,prop=prop,
Imass=x$cacorpo@Imass,Jmass=x$cacorpo@Jmass,a=x$cacorpo@S@Caxes,
b=x$cacorpo@S@Raxes,g=cord2,f=cord1,dmu=diag(x$cacorpo@inertias2[,1]),inertiapc=x$cacorpo@inertias2[,2]),
plottype=plottype,biptype=biptype,pos=pos,arrow=F)
}#end if ellipse
}
