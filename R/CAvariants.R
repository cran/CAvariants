CAvariants <-
function(
  Xtable, mj=NULL,mi=NULL,prop=1,cex=.8, 
 printdims=3, firstaxis=1,lastaxis=2, 
  catype = "CA",plottype="classic",biptype="row",scaleplot=3,posleg="topleft",pos=2,M=min(nrow(Xtable),ncol(Xtable))-1,ell=TRUE ) { 

if (printdims<1) stop(paste("Attention: number of dims for output must be at least 1\n\n"))
if (lastaxis<2) stop(paste("Attention: last axis must be at least 2\n\n"))
if (!any(catype==c("CA","SOCA","DOCA","NSCA","SONSCA","DONSCA"))) stop(paste("Must be CA, DOCA, SOCA, NSCA, SONSCA or DONSCA"))

#if (!any(is.wholenumber(Xtable))) stop(paste("Must be integer values in contingency table"))

# READ DATA FILE
#  assume for now header and row names exist

#Xtable <- read.table(file = datafile, header=header)

#if (header==FALSE) { 
#for (i in 1:dim(Xtable)[1]) rownames(Xtable)[i] <- paste("r",i,sep="")
#for (i in 1:dim(Xtable)[2]) colnames(Xtable)[i] <- paste("c",i,sep="")
#}

X <- as.matrix(Xtable)
#M<-min(nrow(Xtable),ncol(Xtable))-1
  rowlabels <- rownames(Xtable)
  collabels <- colnames(Xtable) 

rows <- dim(X)[1]
cols <- dim(X)[2]
n <- sum(X)
if (is.null(mj)){ 
   mj <- c(1:cols)}
else 
mj<-c(mj) #natural scores for columns  
if (is.null(mi)){ 
   mi <- c(1:rows)}
else 
mi<-c(mi)  #natural scores for rows


maxaxes <- min(rows,cols)-1

S <- switch(catype, "CA"=cabasic(X),  "SOCA"=socabasic(X,mj),"DOCA"=docabasic(X,mi,mj),"NSCA"=nscabasic(X),"SONSCA"=sonscabasic(X,mj),
"DONSCA"=donscabasic(X,mi,mj))

##########################------CA
if(catype=="CA"){
Fmat <- S@RX %*% S@Rweights %*% S@Raxes
Gmat <- S@CX %*% S@Cweights %*% S@Caxes
dmum1 <- diag( (S@mu + (S@mu==0)) * (1-(S@mu==0)) )
Fbi <- S@Cweights %*%  S@Caxes # no orthonormal
Gbi <-S@Rweights %*%   S@Raxes # no orthonormal
pcc <- t(S@CX)
#dimnames(pcc)<-dimnames(X)
inertia <- S@mu*S@mu
comps<-diag(inertia)
Trend<-(Fmat[,firstaxis:lastaxis]%*%t(Gbi[,firstaxis:lastaxis]))

}
#########################################################---------DOCA
if(catype=="DOCA"){
pcc<-S@RX #centered column profile matrix
Gbi <- S@Raxes
Fbi<- S@Caxes 
Gmat <-  S@CX %*%  S@Caxes #row principal coordinates
Fmat <- S@RX  %*% S@Raxes #column principal coordinates
if ((S@Z[1,1]<0) || (S@Z[1,2]<0)){
Gmat<-(-1)*Gmat
Fmat<-(-1)*Fmat
}
reconstruction<-t(Gmat%*%t(S@Cweights%*%Fbi))
dimnames(reconstruction)<-dimnames(X)
inertia <- S@mu #number of inertia of row poly
inertia2<-S@mu2 #number of inertias of column poly
Z<-S@Z
comps <- compstable.exe(Z) 
Icompnames <- c("** Row Components **", "Location", "Dispersion", "Error", "** Column Components **", "Location", "Dispersion", "Error", "** Chi-squared Statistic **")
Jcompnames <- c("Component Value", "P-value")
dimnames(Z) <- list(paste("Poly", 1:(rows - 1)), paste("Poly", 1:(cols - 1)))
dimnames(comps) <- list(paste(Icompnames), paste(Jcompnames))
Trend<-(Fmat[,firstaxis:lastaxis]%*%t(S@Rweights%*%Gbi[,firstaxis:lastaxis]))
#browser()
}
#########################################################################---------SOCA
if(catype=="SOCA"){
pcc<-S@RX
dimnames(pcc)<-dimnames(X)
Gmat <- S@CX %*%  S@Caxes #row principal coordinates
Fmat <- S@RX %*% S@Rweights %*% S@Raxes #column principal coordinates
if (S@Z[1,1]<0){Gmat<-(-1)*Gmat}
Gbi <-S@Raxes
Fbi <- S@Caxes 
inertia <- S@mu*S@mu
inertia2<-S@mu2[-1]
comps<-diag(inertia2)
comps <- compsonetable.exe(S@Z) 
Icompnames <- c( "** Column Components **", "Location", "Dispersion", "Error", "** C-Statistic **")
Jcompnames <- c("Component Value", "P-value")
dimnames(S@Z) <- list(paste("Axis", 1:nrow(S@Z) ), paste("Poly", 1:(cols - 1)))
dimnames(comps) <- list(paste(Icompnames), paste(Jcompnames))
Trend<-(Fmat[,firstaxis:lastaxis]%*%t(S@Rweights%*%Gbi[,firstaxis:lastaxis]))
#browser()
}

####################################-------------NSCA
if(catype=="NSCA"){
Fbi <-  S@Caxes
Gbi <-   S@Raxes 
dmum1 <-diag( (S@mu + (S@mu==0)) * (1-(S@mu==0)) )
pcc<-S@RX
dimnames(pcc)<-dimnames(X)
Gmat <-  S@Raxes %*% dmum1
Fmat <- S@Caxes %*% dmum1
inertia <- S@mu*S@mu
comps<-diag(inertia)
Trend<-(Fmat[,firstaxis:lastaxis]%*%t(S@Rweights%*%Gbi[,firstaxis:lastaxis]))


}
##################################-------------DONSCA
if(catype=="DONSCA"){
Fbi <- S@Caxes
Gbi<-S@Raxes
pcc<-S@RX
dimnames(pcc)<-dimnames(X)
Gmat <-  S@CX %*% S@Cweights %*% S@Caxes #row principal coordinates
Fmat <- S@RX  %*% S@Rweights %*% S@Raxes #column principal coordinates
if (S@Z[1,1]<0){Gmat<-(-1)*Gmat}
inertia <- S@mu
inertia2<-S@mu2 
Z<-sqrt((n-1)*(rows-1))*S@Z
comps <- compstable.exe(Z) 
Icompnames <- c("** Row Components **", "Location", "Dispersion", "Error", "** Column Components **", 
"Location", "Dispersion", "Error", "** Chi-squared Statistic **")
Jcompnames <- c("Component Value", "P-value")
dimnames(Z) <- list(paste("Poly", 1:(rows - 1)), paste("Poly", 1:(cols - 1)))
dimnames(comps) <- list(paste(Icompnames), paste(Jcompnames))
Trend<-(Fmat[,firstaxis:lastaxis]%*%t(Gbi[,firstaxis:lastaxis]))

#browser()
}
############################################------------------SONSCA
if(catype=="SONSCA"){
pcc<-S@RX
dimnames(pcc)<-dimnames(X)
Gmat <- S@CX %*% S@Cweights %*% S@Caxes #column principal coordinates with principal axes
Fmat <- S@RX %*% (S@Rweights) %*% S@Raxes #row principal coordinates with polys
if (S@Z[1,1]<0){Gmat<-(-1)*Gmat}
Gbi <- S@Raxes
Fbi <- S@Caxes 
inertia <- S@mu
inertia2 <- S@mu2
Z<-sqrt((n-1)*(rows-1))*S@Z
comps <- compsonetable.exe(Z) 
Icompnames <- c( "** Column Components **", "Location", "Dispersion", "Error", "** C-Statistic **")
Jcompnames <- c("Component Value", "P-value")
dimnames(S@Z) <- list(paste("Axis", 1:nrow(S@Z)), paste("Poly", 1:(cols - 1)))
dimnames(comps) <- list(paste(Icompnames), paste(Jcompnames))
Trend<-t(Gmat[,firstaxis:lastaxis]%*%t(Fbi[,firstaxis:lastaxis]))

#browser()
}

##################################################################################################################
# OTHER CALCULATIONS

# Calc inertia sum
dmum2 <- diag( 1/(inertia + (S@mu==0)) * (1-(S@mu==0)) )
inertiasum <- sum(inertia)
inertiapc <- 100*inertia/inertiasum
cuminertiapc <- cumsum(inertiapc)
inertiapc <- round(100*inertiapc)/100
cuminertiapc <- round(100*cuminertiapc)/100
inertias <- round(cbind(inertia,inertiapc,cuminertiapc),digits=3)
#browser()

##########################################################
if((catype=="SOCA")|(catype=="SONSCA")|(catype=="DOCA")|(catype=="DONSCA")){
inertiasum2 <- sum(inertia2) #for column categories diag(Z'Z)
inertiapc2 <- 100*inertia2/inertiasum2
cuminertiapc2 <- cumsum(inertiapc2)
inertiapc2 <- round(100*inertiapc2)/100
cuminertiapc2 <- round(100*cuminertiapc2)/100
inertias2 <- round(cbind(inertia2,inertiapc2,cuminertiapc2),digits=3)
}
else inertias2<-inertias
# Calc contributions and correlations

Xstd <- X/sum(X)
if ((catype=="CA")|(catype=="SOCA")|(catype=="DOCA")){
dr <- diag(rowSums(Xstd))}
else {uni<-rep(1,rows)
dr<-diag(uni)}
dc <- diag(colSums(Xstd))
dimnames(Trend)<-list(rowlabels,collabels)

#############################
cacorpo<- new("cacorporateplus",S=S,
DataMatrix=X, rows=rows, cols=cols, 
rowlabels=rowlabels, collabels=collabels,
Rprinccoord=Fmat, Cprinccoord=Gmat, Rstdcoord=Fbi, Cstdcoord=Gbi,
 inertiasum=inertiasum, inertias=inertias, inertias2=inertias2,comps=comps,
  printdims=printdims,maxaxes=maxaxes,catype=catype,mj=mj,mi=mi,pcc=pcc,Jmass=dc,Imass=dr,Trend=Trend)

#browser()
printcacorporateplus(cacorpo)

plotcacorporateplus(cacorpo,cex=cex,firstaxis=firstaxis,lastaxis=lastaxis,inert=inertias,inertsum=inertiasum,prop=prop,M=M,catype=catype,
biptype=biptype,plottype=plottype,scaleplot=scaleplot,posleg=posleg,pos=pos,ell=ell)

}
