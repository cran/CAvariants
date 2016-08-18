compsonetable.exe<- 
function(Z){
nr<-nrow(Z)
nc<-ncol(Z)
tZZ <- t(Z) %*% (Z)
factor <- sum(diag(tZZ))
comps <- matrix(0, nrow=3, ncol=2)
###############################
#     Column Category         #
###############################
for (j in 1:3){
comps[j, 1] <- tZZ[j, j]#Location Component for Category 2#
comps[j, 2] <- 1 - pchisq(comps[j, 1], ncol(Z))#P-value of the Location Comp for Category 2#
}
compsClast1 <- factor -sum (comps[, 1])#Error of Row Components for Category 1#
if(ncol(Z) > 3) {
compsClast2 <- 1 - pchisq(compsClast1, (nc-3) * nr )#P-value for the Errors of Category 1#
Ccol=rbind(comps,c(compsClast1,compsClast2))
}
else {
compsClast2 <- 0
Ccol=rbind(comps,c(compsClast1,compsClast2))
}
compsCtot1 <- factor
compsCtot2 <- 1 - pchisq(compsCtot1, nr * nc)
Ccol=rbind(Ccol,c(compsCtot1,compsCtot2))
list(comps=as.matrix(Ccol))
}
