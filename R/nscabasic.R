nscabasic <-
function(Xtable) { 
X <- Xtable/sum(Xtable)
#r<- min(dim(X))-1
I<-nrow(X)
J<-ncol(X)
Imass<-rowSums(X)
tauden <- 1 - sum(Imass^2)
rsums <- as.vector(rowSums(X))
csums <- as.vector(colSums(X))
drm1 <- diag( 1/( rsums + (rsums==0) ) * (1-(rsums==0)) )
dcm1 <- diag( 1/( csums + (csums==0) ) * (1-(csums==0)) )
drmh<-diag(rep(1,I)) #change the metric in NSCA
dcmh <- sqrt(dcm1)
#Z <- 1/sqrt(tauden)*(drmh %*% ( X - rsums %*% t(csums) ) %*% dcmh)#tau index
Z <- (drmh %*% ( X - rsums %*% t(csums) ) %*% dcmh) #only numerator
#Yeigu<-eigen(Z%*%t(Z))
#Caxes<-Yeigu$vectors
#Yeigv<-eigen(t(Z)%*%Z)
#Raxes<-dcmh%*%Yeigv$vectors
Y <- svd(Z,nu=I,nv=J)
Raxes<-dcmh%*%Y$v
Caxes<-Y$u
mu <- Y$d
#tau<-sum(mu^2)/tauden
R <- drm1 %*% X
C <- dcm1 %*% t(X)
#browser()
#NSCA<-new("cabasicresults",
 #         RX=R,CX=C,Rweights=dcmh,Cweights=drmh,
  #        Raxes=dcmh%*%Y$v,Caxes=Y$u,mu=mu,tauDen=tauden,catype="NSCA")
#----------------------------
#browser()
resnsca<-(list(
          RX=R,CX=C,Rweights=dcmh,Cweights=drmh,
          Raxes=Raxes,Caxes=Caxes,mu=mu,tauDen=tauden,catype="NSCA"))
return(resnsca)
}
