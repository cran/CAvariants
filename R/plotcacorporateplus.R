plotcacorporateplus <-
function(x,  firstaxis=1, lastaxis=2, cex=1,prop=1,inert,inertsum,Xdata,catype="CA",
biptype="rows",plottype="classic",M=(min(nrow(x),ncol(x))-1),scaleplot=3,posleg="topleft",pos=2,ell=TRUE) {

## internal function to plot a  single picture
##

if ((firstaxis<1)|(firstaxis>x@maxaxes-1)) stop(paste("incorrect first axis =", firstaxis, "\n\n")) 
if (lastaxis>x@maxaxes) stop(paste("incorrect last axis =", lastaxis, "\n\n")) 
if (firstaxis>=lastaxis) stop(paste("last axis must be greater than first axis\n\n"))

tworowS <- rowSums(x@DataMatrix>0)==2
twocolS <- colSums(x@DataMatrix>0)==2

# Groups file must have no blank line at start and only one between sections
# group number   group name   symbol   colour   plot ellipse? 

rowgroup <- list(1:x@rows,rep(1,x@rows))
rowgrlab <- list(1,"","*","red","T")
colgroup <- list(1:x@cols,rep(1,x@cols))
colgrlab <- list(1,"","+","blue","T")
######################################################
# Plot row and col coordinates
#########################################
if ((plottype=="Classical")|(plottype=="classical")|(plottype=="classic")|(plottype=="c")) {
nthings<-x@cols
nvars<-x@rows
cord1<- x@Cprinccoord 
cord2<-x@Rprinccoord
inertiapc=x@inertias[,2] #inertia in percentage of row axes 
dimnames(cord1)[1]<-dimnames(x@DataMatrix)[2]
dimnames(cord2)[1]<-dimnames(x@DataMatrix)[1]
thinggroup<-colgroup
thinggrlab<-colgrlab
vargroup<-rowgroup
vargrlab<-rowgrlab
thinglabels<-x@collabels
varlabels<-x@rowlabels
main="Classical plot"
if ((x@catype=="DONSCA")|(x@catype=="DOCA")|(x@catype=="SOCA")|(x@catype=="SONSCA"))
{
cat("\n NO CLASSICAL PLOT ONLY BIPLOT (coordinates are biplot coordinates, please change plottype)\n")
cord1<-x@Cstdcoord/scaleplot
cord2<-x@Rprinccoord*scaleplot
inertiapc=x@inertias2[,2] #inertia in percentage of row axes 
}

}

         if ((plottype=="Biplot")|(plottype=="biplot")|(plottype=="bip")|(plottype=="b")){
if ((biptype=="rows")|(biptype=="Rows")|(biptype=="row")|(biptype=="r")) 
{
plottype<-"biplot"
#biptype<-"row"
cord1<-x@Rprinccoord*scaleplot
cord2<-x@Cstdcoord/scaleplot
nthings<-x@rows
nvars<-x@cols
thinglabels<-x@rowlabels
varlabels<-x@collabels
thinggroup<-rowgroup
thinggrlab<-rowgrlab
vargroup<-colgroup
vargrlab<-colgrlab
main<-"Row Isometric Biplot"
inertiapc=x@inertias[,2] #inertia of column poly
dimnames(cord2)[1]<-dimnames(x@DataMatrix)[2]
dimnames(cord1)[1]<-dimnames(x@DataMatrix)[1]

if ((x@catype=="DONSCA")|(x@catype=="DOCA")|(x@catype=="SOCA")|(x@catype=="SONSCA"))
{
cord2<-x@Rprinccoord*scaleplot
cord1<-x@Cstdcoord/scaleplot
nthings<-x@cols
nvars<-x@rows
thinglabels<-x@collabels
varlabels<-x@rowlabels
thinggroup<-colgroup
thinggrlab<-colgrlab
vargroup<-rowgroup
vargrlab<-rowgrlab
inertiapc=x@inertias2[,2] #inertia of column poly
dimnames(cord2)[1]<-dimnames(x@DataMatrix)[1]
dimnames(cord1)[1]<-dimnames(x@DataMatrix)[2]

}#end catype

} #end bip row
else{
plottype<-"biplot"
#biptype<-"columns"
cord1<- x@Cprinccoord*scaleplot
cord2<-x@Rstdcoord/scaleplot
nthings<-x@cols
nvars<-x@rows
thinggroup<-colgroup
thinggrlab<-colgrlab
vargroup<-colgroup
vargrlab<-rowgrlab
thinglabels<-x@collabels
varlabels<-x@rowlabels
main<-"Column Isometric Biplot"
inertiapc=x@inertias[,2] #inertia of row 
dimnames(cord1)[1]<-dimnames(x@DataMatrix)[2]
dimnames(cord2)[1]<-dimnames(x@DataMatrix)[1]

if ((x@catype=="DONSCA")|(x@catype=="DOCA")|(x@catype=="SOCA")|(x@catype=="SONSCA"))
{
cord2<- x@Cprinccoord*scaleplot
cord1<-x@Rstdcoord/scaleplot
nthings<-x@rows
nvars<-x@cols
thinggroup<-rowgroup
thinggrlab<-rowgrlab
vargroup<-colgroup
vargrlab<-colgrlab
thinglabels<-x@rowlabels
varlabels<-x@collabels

inertiapc=x@inertias[,2] #inertia of row poly
dimnames(cord1)[1]<-dimnames(x@DataMatrix)[1]
dimnames(cord2)[1]<-dimnames(x@DataMatrix)[2]

}#end catype

}#end bip column
}


###################################################################################ok without choice plottype
#repeat{

if ((catype=="DOCA")|(catype=="SOCA")|(catype=="SONSCA")|(catype=="DONSCA"))
{
#dimnames(x@pcc)<-list(x@rowlabels,x@collabels)
 cat("\n Looking at the Trends of rows and columns\n")
#browser()
#trendplot(x@mj,(x@pcc), position="topleft",main="Original Trends of Centered #Column Profile",xlab="ordered scores")
#windows()
#trendplot(x@mi,t(x@pcc), position="topleft",main="Original Trends of Centered #Column Profile",xlab="ordered scores")
############################################################## reconstructed TREND
#windows()
trendplot(x@mj,(x@Trend), posleg=posleg,main="Reconstructed rows of the centred column profile",xlab="ordered scores",prop=prop)
dev.new()
trendplot(x@mi,t(x@Trend), posleg=posleg,main="Reconstructed columns of the centred column profile",xlab="ordered scores",prop=prop)
#browser()
}

##############################################
picsize1<-c(range(cord1[,c(firstaxis,lastaxis)], cord2[,c(firstaxis,lastaxis)])/prop)
picsizeR<-c(range(cord1[,c(firstaxis,lastaxis)], cord1[,c(firstaxis,lastaxis)])/prop)
picsize2<-c(range(cord2[,c(firstaxis,lastaxis)], cord2[,c(firstaxis,lastaxis)])/prop)

if (picsize1[1]>=picsize1[2]) stop(paste("incorrect axis scale picsize =", picsize1[1], picsize1[2], "\n\n"))
########################################################################################## 
if ((x@catype=="DONSCA")||(x@catype=="DOCA"))
{
plotone (firstaxis,lastaxis,plottype=plottype,things=x@catype,nthings,nvars,cord1,cord2,
inertiapc=inertiapc, thinggroup,thinggrlab,vargroup,vargrlab, thinglabels, varlabels,picsize=picsize1,cex=cex,type="b",catype="DOCA",pos=pos) 
}
######################
if ((x@catype=="SOCA")||(x@catype=="SONSCA"))
{
if (biptype=="row"){type="b"}
else {type="p"}
plotone (firstaxis,lastaxis,plottype=plottype,things=x@catype,nthings,nvars,cord1,cord2,
inertiapc=inertiapc, thinggroup,thinggrlab,vargroup,vargrlab, thinglabels, varlabels,picsize=picsize1,cex=cex,catype="SOCA",type=type,pos=pos) 
}

###############################
if ((x@catype=="CA")||(x@catype=="NSCA"))
{
plotone (firstaxis,lastaxis,plottype=plottype,things=x@catype,nthings,nvars,cord1,cord2,
inertiapc=inertiapc, thinggroup,thinggrlab,vargroup,vargrlab, thinglabels, varlabels,picsize=picsize1,cex=cex,type="p",catype="CA",pos=pos) 
}
################################################################################
#cat("\nIncluding Beh's Confidence Ellipses\n")
#####################################################
if (ell==TRUE) {
if ((catype=="CA")&(plottype=="biplot")&(biptype=="row")){
cordr<-cord2
cordc<-cord1
cord1<-cordr
cord2<-cordc
}
if ((catype=="NSCA")&(plottype=="biplot")&(biptype=="row")){
cordr<-cord2
cordc<-cord1
cord1<-cordr
cord2<-cordc
}

if (((catype=="DOCA")|(catype=="SOCA")|(catype=="SONSCA")|(catype=="DONSCA")) & (plottype=="biplot")&((biptype=="column")|(biptype=="col"))){
cordr<-cord2
cordc<-cord1
cord1<-cordr
cord2<-cordc
}
dev.new()
M<-(min(nrow(x@DataMatrix), ncol(x@DataMatrix)) -1)

switch(x@catype, "CA"=caellipse(N=x@DataMatrix,a1=firstaxis,a2=lastaxis,M=M,prop=prop,
Imass=x@Imass,Jmass=x@Jmass,a=x@S@Caxes,b=x@S@Raxes,g=cord1,f=cord2,dmu=diag(x@inertias[,1]),inertiapc=x@inertias[,2],
plottype=plottype,biptype=biptype,pos=pos,arrow=T), 

"SOCA"=caellipse(N=x@DataMatrix,a1=firstaxis,a2=lastaxis,M=M,prop=prop,
Imass=x@Imass,Jmass=x@Jmass,a=x@S@Caxes,b=x@S@Raxes,g=cord1,f=cord2,dmu=diag(x@inertias2[,1]),inertiapc=x@inertias2[,2],
plottype=plottype,biptype=biptype,pos=pos,arrow=F), 

"DOCA"=caellipse(N=x@DataMatrix,a1=firstaxis,a2=lastaxis,M=(M-1),prop=prop,
Imass=x@Imass,Jmass=x@Jmass,a=x@S@Caxes,b=x@S@Raxes,f=cord2,g=cord1,dmu=diag(x@inertias2[,1]),inertiapc=x@inertias2[,2],
plottype=plottype,biptype=biptype,pos=pos,arrow=F), 

"NSCA"=nsca.ellipse(x@DataMatrix,a1=firstaxis,a2=lastaxis,M=M,prop=prop,
Imass=x@Imass,Jmass=x@Jmass,a=x@S@Caxes,b=x@S@Raxes,g=cord1,f=cord2,dmu=diag(x@inertias[,1]),inertiapc=x@inertias[,2],
plottype=plottype,biptype=biptype,pos=pos,arrow=T),

"SONSCA"=nsca.ellipse(x@DataMatrix,a1=firstaxis,a2=lastaxis,M= M,prop=prop,
Imass=x@Imass,Jmass=x@Jmass,a=x@S@Caxes,b=x@S@Raxes,f=cord2,g=cord1,dmu=diag(x@inertias2[,1]),inertiapc=x@inertias2[,2],
plottype=plottype,biptype=biptype,pos=pos,arrow=F), 

"DONSCA"=nsca.ellipse(x@DataMatrix,a1=firstaxis,a2=lastaxis,M= (M-1),prop=prop,
Imass=x@Imass,Jmass=x@Jmass,a=x@S@Caxes,b=x@S@Raxes,f=cord2,g=cord1,dmu=diag(x@inertias2[,1]),inertiapc=x@inertias2[,2]),
plottype=plottype,biptype=biptype,pos=pos,arrow=F)
}#end if ellipse
}#end plotfunction
