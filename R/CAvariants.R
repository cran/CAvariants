CAvariants<-
function(
  Xtable, mj=NULL,mi=NULL,firstaxis=1,lastaxis=2,
  catype = "CA",M=min(nrow(Xtable),ncol(Xtable))-1,alpha=0.05) { 
#if (printdims<1) stop(paste("Attention: number of dims for output must be at least 1\n\n"))
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
r<- min(rows,cols)-1
S <- switch(catype, "CA"=cabasic(X),  "SOCA"=socabasic(X,mj),"DOCA"=docabasic(X,mi,mj),"NSCA"=nscabasic(X),"SONSCA"=sonscabasic(X,mj),
"DONSCA"=donscabasic(X,mi,mj))
##########################------CA
if(catype=="CA"){
Fmat <- S$RX %*% S$Rweights %*% S$Raxes[,1:r]
Gmat <- S$CX %*% S$Cweights %*% S$Caxes[,1:r]
#dmum1 <- diag( (S$mu + (S$mu==0)) * (1-(S$mu==0)) )
Fbi <- S$Cweights %*%  S$Caxes[,1:r] # no orthonormal
Gbi <-S$Rweights %*%   S$Raxes[,1:r] # no orthonormal
pcc <- t(S$CX)
#dimnames(pcc)<-dimnames(X)
tau=NULL
tauden=NULL
inertia <- (S$mu[1:r]^2) #please check!!
inertiasum <- sum(S$mu^2)
inertiasum2 <- sum(S$mu^2)
t.inertia<-inertiasum*n
comps<-diag(inertia)
Trend<-(Fmat[,firstaxis:lastaxis]%*%t(Gbi[,firstaxis:lastaxis]))
Z<-Trend
dimnames(Z) <- list(rowlabels,collabels)
}
#########################################################---------DOCA
if(catype=="DOCA"){
Z<-S$Z/sqrt(n)
pcc<-S$RX #centered column profile matrix
Gbi <- S$Raxes[,1:r]
Fbi<- S$Caxes[,1:r] 
Gbi2 <- S$Raxes
Fbi2<- S$Caxes
Gmat <-  S$CX %*%  S$Caxes[,1:r] #row principal coordinates
Fmat <- S$RX  %*% S$Raxes[,1:r] #column principal coordinates
nr<-nrow(Z)
nc<-ncol(Z)
#if (nr>2){
#if  ((Z[2,2]<0) & (Z[1,2]>0)||(Z[2,2]>0) & (Z[1,2]<0)){
#Gmat<-(-1)*Gmat
#Fmat<-(-1)*Fmat
#}}
reconstruction<-t(Gmat%*%t(S$Cweights%*%Fbi))
dimnames(reconstruction)<-dimnames(X)
inertia <- S$mu[1:r]/n #number of inertia of row poly
inertia2<-S$mu2[1:r]/n #number of inertias of column poly
inertiasum <- sum(S$mu)/n
inertiasum2 <- sum(S$mu2)/n
t.inertia<-inertiasum*n
Z1<-S$Z
tau=NULL
tauden=NULL
comps <- compstable.exe(Z1)
#----------------------------------------rows comps
if (nr>3){
Rcompnames <- c( "Location", "Dispersion", "Cubic","Error", "** Chi-squared Statistic **")}
if(nr == 3){
Rcompnames <- c( "Location", "Dispersion","Error","** Chi-squared Statistic **")}
if(nr == 2){
Rcompnames <- c( "Location", "Dispersion","** Chi-squared Statistic **")}

#----------------------columns comps
if (nc>3){
Ccompnames <- c( "Location", "Dispersion", "Cubic","Error", "** Chi-squared Statistic **")}
if (nc==3){
Ccompnames <- c( "Location", "Dispersion","Error","** Chi-squared Statistic **")}
if (nc==2){
Ccompnames <- c( "Location", "Dispersion","** Chi-squared Statistic **")}
Jcompnames <- c("Inertia Value", "P-value")
dimnames(Z) <- list(paste("u", 1:nr,sep=""), paste("v", 1:nc,sep=""))
#browser()
dimnames(comps$compsR) <- list(paste(Rcompnames), paste(Jcompnames))
dimnames(comps$compsC) <- list(paste(Ccompnames), paste(Jcompnames))
Trend<-(Fmat[,firstaxis:lastaxis]%*%t(S$Rweights%*%Gbi[,firstaxis:lastaxis]))
#browser()
}
#########################################################################---------SOCA
if(catype=="SOCA"){
pcc<-S$RX
Z<-S$Z/sqrt(n)
nr<-nrow(Z)
nc<-ncol(Z)
dimnames(pcc)<-dimnames(X)
Gmat <- S$CX %*%  S$Caxes[,1:r] #row principal coordinates
Fmat <- S$RX %*% S$Rweights %*% S$Raxes[,1:r] #column principal coordinates
#if (nr>2){
#if  ((Z[2,2]<0) & (Z[1,2]>0)||(Z[2,2]>0) & (Z[1,2]<0)){
#Gmat<-(-1)*Gmat
#Fmat<-(-1)*Fmat
#}}
Gbi <-S$Rweights %*%S$Raxes[,1:r]
Fbi <- S$Cweights %*%S$Caxes[,1:r] 
inertia <- (S$mu[1:r]^2)/n
inertia2<-(S$mu2[1:r])/n
inertiasum <- sum(S$mu^2)/n
inertiasum2 <- sum(S$mu2)/n
t.inertia<-inertiasum*n
tauden=NULL
tau=NULL
Z1<-S$Z
comps <- compsonetable.exe(Z1) 
#----------------------columns comps
if (nc>3){
Ccompnames <- c( "Location", "Dispersion", "Cubic","Error", "** Chi-squared Statistic **")}
if (nc==3){
Ccompnames <- c( "Location", "Dispersion","Error","** Chi-squared Statistic **")}
if (nc==2){
Ccompnames <- c( "Location", "Dispersion","** Chi-squared Statistic **")}
Jcompnames <- c("Inertia Value", "P-value")
dimnames(Z) <- list(paste("m", 1:nr,sep="" ), paste("v", 1:nc,sep=""))
#browser()
dimnames(comps$comps) <- list(Ccompnames, Jcompnames)
Trend<-(Fmat[,firstaxis:lastaxis]%*%t(S$Rweights%*%Gbi[,firstaxis:lastaxis]))
#-----------------------------
}
####################################-------------NSCA
if(catype=="NSCA"){
Fbi <-  S$Caxes[,1:r]
Gbi <-   S$Raxes[,1:r] 
#dmum1 <-diag( (S$mu + (S$mu==0)) * (1-(S$mu==0)) )
dmum1 <-diag( S$mu [1:r])
pcc<-S$RX
dimnames(pcc)<-dimnames(X)
Gmat <-  S$Raxes[,1:r] %*% dmum1
Fmat <- S$Caxes[,1:r] %*% dmum1
tauden<-S$tauDen
inertia <- S$mu[1:r]^2
inertiasum <- sum(S$mu^2)
inertiasum2 <- sum(S$mu^2)
tau<-sum(inertia)/tauden
 t.inertia <- (sum(S$mu^2)/tauden)*(n-1)*(rows-1)
comps<-diag(inertia)
Trend<-(Fmat[,firstaxis:lastaxis]%*%t(S$Rweights%*%Gbi[,firstaxis:lastaxis]))
Z<-Trend
dimnames(Z) <- list(rowlabels,collabels)
}
##################################-------------DONSCA
if(catype=="DONSCA"){
Fbi <- S$Caxes[,1:r]
Gbi<-S$Raxes[,1:r]
pcc<-S$RX
Z<-S$Z
nr<-nrow(Z)
nc<-ncol(Z)
dimnames(pcc)<-dimnames(X)
Gmat <-  S$CX %*% S$Cweights %*% S$Caxes[,1:r] #row principal coordinates
Fmat <- S$RX  %*% S$Rweights %*% S$Raxes[,1:r] #column principal coordinates
#if (nr>2){
#if  ((Z[2,2]<0) & (Z[1,2]>0)||(Z[2,2]>0) & (Z[1,2]<0)){
#Gmat<-(-1)*Gmat
#Fmat<-(-1)*Fmat
#}}
inertia <- S$mu[1:r]
inertia2<-S$mu2[1:r] 
inertiasum <- sum(S$mu)
inertiasum2 <- sum(S$mu2)
tauden<-S$tauDen
Z2<-1/sqrt(tauden)*sqrt((n-1)*(rows-1))*S$Z
#Z2<-sqrt((n-1)*(rows-1))*S$Z #when tau
 t.inertia <- (sum(S$mu)/tauden)*(n-1)*(rows-1)
tau<-sum(inertia)/tauden
comps <- compstable.exe(Z2) 
#----------------------------------------rows comps
if (nr>3){
Rcompnames <- c( "Location", "Dispersion", "Cubic","Error", "** Chi-squared Statistic **")}
if(nr == 3){
Rcompnames <- c( "Location", "Dispersion","Error","** Chi-squared Statistic **")}
if(nr == 2){
Rcompnames <- c( "Location", "Dispersion","** Chi-squared Statistic **")}

#----------------------columns comps
if (nc>3){
Ccompnames <- c( "Location", "Dispersion", "Cubic","Error", "** Chi-squared Statistic **")}
if (nc==3){
Ccompnames <- c( "Location", "Dispersion","Error","** Chi-squared Statistic **")}
if (nc==2){
Ccompnames <- c( "Location", "Dispersion","** Chi-squared Statistic **")}
Jcompnames <- c("Inertia Value", "P-value")
dimnames(Z) <- list(paste("u", 1:(rows-1 ),sep=""), paste("v", 1:(cols -1),sep=""))
#dimnames(Z) <- list(paste("u", 1:(rows ),sep=""), paste("v", 1:(cols -1),sep=""))
dimnames(comps$compsR) <- list(paste(Rcompnames), paste(Jcompnames))
dimnames(comps$compsC) <- list(paste(Ccompnames), paste(Jcompnames))
Trend<-(Fmat[,firstaxis:lastaxis]%*%t(Gbi[,firstaxis:lastaxis]))
}
############################################------------------SONSCA
if(catype=="SONSCA"){
pcc<-S$RX
dimnames(pcc)<-dimnames(X)
Z<-S$Z
nr<-nrow(Z)
nc<-ncol(Z)
Gmat <- S$CX %*% S$Cweights %*% S$Caxes[,1:r] #column principal coordinates with principal axes
Fmat <- S$RX %*% (S$Rweights) %*% S$Raxes[,1:r] #row principal coordinates with polys
Gbi <- S$Raxes[,1:r]
Fbi <- S$Caxes[,1:r] 
inertia <- S$mu[1:r]
inertia2 <- S$mu2[1:r]
inertiasum <- sum(S$mu) #computed through diag(ZZ')
inertiasum2 <- sum(S$mu2) #for column categories diag(Z'Z)
tauden<-S$tauDen
tau<-sum(inertia)/tauden
Z1<-1/sqrt(tauden)*sqrt((n-1)*(rows-1))*S$Z
 t.inertia <- (inertiasum/tauden)*(n-1)*(rows-1)
comps <- compsonetable.exe(Z1) 
#----------------------columns comps
if (nc>3){
Ccompnames <- c( "Location", "Dispersion", "Cubic","Error", "** Chi-squared Statistic **")}
if (nc==3){
Ccompnames <- c( "Location", "Dispersion","Error","** Chi-squared Statistic **")}
if (nc==2){
Ccompnames <- c( "Location", "Dispersion","** Chi-squared Statistic **")}
Jcompnames <- c("Inertia Value", "P-value")
dimnames(Z) <- list(paste("m", 1:nr,sep=""), paste("v", 1:nc,sep=""))
dimnames(comps$comps) <- list(Ccompnames, Jcompnames)
Trend<-t(Gmat[,firstaxis:lastaxis]%*%t(Fbi[,firstaxis:lastaxis]))
}

##################################################################################################################
# OTHER CALCULATIONS

# Calc inertia sum
inertiapc <- 100*inertia/inertiasum
cuminertiapc <- cumsum(inertiapc)
inertiapc <- (100*inertiapc)/100
cuminertiapc <- (100*cuminertiapc)/100
inertias <- cbind(inertia,inertiapc,cuminertiapc)

##########################################################
if((catype=="SOCA")|(catype=="SONSCA")|(catype=="DOCA")|(catype=="DONSCA")){
#inertiasum2 <- sum(inertia2) #for column categories diag(Z'Z)
inertiapc2 <- 100*inertia2/inertiasum2
cuminertiapc2 <- cumsum(inertiapc2)
inertiapc2 <- (100*inertiapc2)/100
cuminertiapc2 <- (100*cuminertiapc2)/100
inertias2 <- cbind(inertia2,inertiapc2,cuminertiapc2)
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
dimnames(dr)<-list(rowlabels,rowlabels)
dimnames(dc)<-list(collabels,collabels)
dimnames(inertias)[[1]]<-paste("value", 1:r,sep="")
dimnames(inertias2)[[1]]<-paste("value", 1:r,sep="")

dimnames(Fmat)<-list(rowlabels,paste("axes", 1:r,sep=""))
dimnames(Fbi)<-list(rowlabels,paste("axes", 1:r,sep=""))
dimnames(Gmat)<-list(collabels,paste("axes", 1:r,sep=""))
dimnames(Gbi)<-list(collabels,paste("axes", 1:r,sep=""))
names(mi)<-rowlabels
names(mj)<-collabels
#---------------------------------------------------
cord1<-Gmat #check here
cord2<-Fmat
if ((catype=="DOCA")|(catype=="SOCA")|(catype=="SONSCA")|(catype=="DONSCA")){
cordr<-cord2
cordc<-cord1
cord1<-cordr
cord2<-cordc
}
#-------------------------------------------------------------
#for printing the characteristics of  ellipses 
#----------------------------------------------------
#if ((ellcomp==TRUE)&&(catype=="CA")|(ellcomp==TRUE)&&(catype=="NSCA")){
#if (ellcomp==TRUE){
#---------------------------------p-values
inertiapc<-inertias[,2]
cord1<-Fmat
cord2<-Gmat
a<-Fbi
b<-Gbi
 Inames<-rowlabels
Jnames<-collabels
I<-rows
J<-cols
dmu<-sqrt(inertias[,1])
#t.inertia<-inertiasum*n
Imass<-dr
Jmass<-dc
#if (catype=="NSCA"){
#    t.inertia <- (sum(dmu^2)/tauden)*(n-1)*(I-1)
#  }
#-----------------------------axis ellipses
chisq.val <- qchisq(1 - alpha, df = (I - 1) * (J - 1))
hlax1.row <- vector(mode = "numeric", length = I)
hlax2.row <- vector(mode = "numeric", length = I)
hlax1.col <- vector(mode = "numeric", length = J)
hlax2.col <- vector(mode = "numeric", length = J)
#browser()
if (M > 2) {
  for (i in 1:I) {
    hlax1.row[i] <- dmu[1] * sqrt(abs((chisq.val/(t.inertia)) * 
                                          (1/Imass[i,i] - sum(a[i, 3:M]^2))))
    hlax2.row[i] <- dmu[2] * sqrt(abs((chisq.val/(t.inertia)) * 
                                          (1/Imass[i,i] - sum(a[i, 3:M]^2))))
  }
  for (j in 1:J) {
    hlax1.col[j] <- dmu[1] * sqrt(abs((chisq.val/(t.inertia)) * 
                                          (1/Jmass[j,j] - sum(b[j, 3:M]^2))))
    hlax2.col[j] <- dmu[2] * sqrt(abs((chisq.val/(t.inertia)) * 
                                          (1/Jmass[j,j] - sum(b[j, 3:M]^2))))
      }
            }
if (M == 2) {
    for (i in 1:I) {
    hlax1.row[i] <- dmu[1] * sqrt(abs((chisq.val/(t.inertia)) * 
                                          (1/(Imass)[i,i])))
    hlax2.row[i] <- dmu[2] * sqrt(abs((chisq.val/(t.inertia)) * 
                                          (1/(Imass)[i,i])))
  }
  for (j in 1:J) {
    hlax1.col[j] <- dmu[1] * sqrt(abs((chisq.val/(t.inertia)) * 
                                          (1/(Jmass)[j,j])))
    hlax2.col[j] <- dmu[2] * sqrt(abs((chisq.val/(t.inertia)) * 
                                          (1/(Jmass)[j,j])))
  }
}
eccentricity <- abs(1 - (dmu[2]^2/dmu[1]^2))^(1/2)
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
    pvalrow[i] <- 1- pchisq(t.inertia*((1/Imass[i,i] - sum(a[i, 
                                                                 3:M]^2))^(-1)) * (cord1[i, 1]^2/dmu[1]^2 + cord1[i, 
                                                                                                               2]^2/dmu[2]^2), df = (I - 1) * (J - 1))
  }
  else {
    pvalrow[i] <- 1-pchisq(t.inertia* (Imass[i,i]) * 
                             (cord1[i, 1]^2/dmu[1]^2 + cord1[i, 2]^2/dmu[2]^2), 
                           df = (I - 1) * (J - 1))
  }
}
for (j in 1:J) {
  if (M > 2) {
    pvalcol[j] <-  1-pchisq(t.inertia* ((1/Jmass[j,j] - 
                                              sum(b[j, 3:M]^2))^(-1)) * (cord2[j, 1]^2/dmu[1]^2 + 
                                                                           cord2[j, 2]^2/dmu[2]^2), df = (I - 1) * (J - 1))
  }
  else {
    pvalcol[j] <- 1- pchisq(t.inertia * (Jmass[j,j]) * 
                              (cord2[j, 1]^2/dmu[1]^2 + cord2[j, 2]^2/dmu[2]^2), 
                            df = (I - 1) * (J - 1))
  }
}
summ.name <- c("HL Axis 1", "HL Axis 2", "Area", "P-value")
row.summ <- cbind(hlax1.row, hlax2.row, area.row, pvalrow)
dimnames(row.summ) <- list(paste(Inames), paste(summ.name))
col.summ <- cbind(hlax1.col, hlax2.col, area.col, pvalcol)
dimnames(col.summ) <- list(paste(Jnames), paste(summ.name))
#}
#if ((ellcomp==FALSE)|(ellcomp==TRUE)&&(catype=="DOCA")|(ellcomp==TRUE)&&(catype=="DONSCA")|(ellcomp==TRUE)&&(catype=="SOCA")|(ellcomp==TRUE)&&(catype=="SONSCA")){
#if (ellcomp==FALSE){
#eccentricity=NULL
#row.summ=NULL
#col.summ=NULL
#}
#--------------------------------------------------
resultCA<-list(Xtable=X, rows=rows, cols=cols, r=r,n=n,
rowlabels=rowlabels, collabels=collabels,
Rprinccoord=Fmat, Cprinccoord=Gmat, Rstdcoord=Fbi, Cstdcoord=Gbi,tauden=tauden,tau=tau,
inertiasum2=inertiasum2, inertiasum=inertiasum, inertias=inertias, inertias2=inertias2,t.inertia=t.inertia,comps=comps,
 catype=catype,mj=mj,mi=mi,pcc=pcc,Jmass=dc,Imass=dr,
Innprod=Trend,Z=Z,M=M,eccentricity=eccentricity,row.summ=row.summ,col.summ=col.summ)
class(resultCA)<-"CAvariants"
return(resultCA)
#-----------------------------------------------
}
