nscabasic <-
function(X) { 
X <- X/sum(X)
r<- min(dim(X))-1
Imass<-rowSums(X)
tauden <- 1 - sum(Imass^2)

rsums <- as.vector(rowSums(X))
csums <- as.vector(colSums(X))
drm1 <- diag( 1/( rsums + (rsums==0) ) * (1-(rsums==0)) )
dcm1 <- diag( 1/( csums + (csums==0) ) * (1-(csums==0)) )
drmh<-diag(rep(1,nrow(X))) #change the metric in NSCA
dcmh <- sqrt(dcm1)
Z <- 1/sqrt(tauden)*(drmh %*% ( X - rsums %*% t(csums) ) %*% dcmh)
Y <- svd(Z)
mu <- Y$d
tau<-sum(mu^2)
R <- drm1 %*% X
C <- dcm1 %*% t(X)
#browser()
NSCA<-new("cabasicresults",
          RX=R,CX=C,Rweights=dcmh,Cweights=drmh,
          Raxes=dcmh%*%Y$v,Caxes=Y$u,r=r,mu=mu,tau=tau,tauDen=tauden,catype="NSCA")
}
