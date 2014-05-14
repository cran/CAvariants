socabasic <-
function (xo,mj) 
{
n<-sum(xo)
    x <- xo/n
rsums <- as.vector(rowSums(x))
csums <- as.vector(colSums(x))
di<-diag(rsums)
dj<-diag(csums)
    Bpoly <- emerson.poly(mj, csums)
#  Bpoly <- orthopoly.exe(c(csums))[,-1]

    Bpoly2 <- sqrt(dj) %*% Bpoly
######################################################
drm1 <- diag( 1/( rsums + (rsums==0) ) * (1-(rsums==0)) )
dcm1 <- diag( 1/( csums + (csums==0) ) * (1-(csums==0)) )
drmh <- sqrt(drm1)
dcmh <- sqrt(dcm1)
ratio <- drmh %*% ( x - rsums %*% t(csums) ) %*% dcmh*sqrt(n)
u<-svd(ratio)$u
mu <- svd(ratio)$d
R <- drm1 %*% x
C <- dcm1 %*% t(x)
Z <- t(u) %*% ratio %*% Bpoly2 #useful to check coordinates
ZtZ <- Z%*%t(Z)
tZZ <- t(Z)%*%Z
r <- sum(mu>1e-15)
rmax <- min(dim(xo)) - 1
    r <- sum(mu > 1e-15)
    if (r < rmax) {
        mu[(r + 1):rmax] <- 0
        Raxes[, (r + 1):rmax] <- 0
        Caxes[, (r + 1):rmax] <- 0
    }
soca <- new("cabasicresults",
          RX=R,CX=C,Rweights=dcmh,Cweights=drmh,
         Raxes= Bpoly,Caxes=u,r=r,mu=mu,mu2=diag(ZtZ),tau=0,tauDen=0,catype="SOCA",Z=Z,ZtZ=ZtZ,tZZ=tZZ)

}
