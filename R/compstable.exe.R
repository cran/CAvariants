compstable.exe<-
function(Z){
tZZ <- t(Z) %*% Z
ZtZ <- Z %*% t(Z)
factor <- sum(diag(tZZ))
nr<-nrow(Z)
nc<-ncol(Z)
compsR <- matrix(0, nrow = 3, ncol = 2)
compsC <- matrix(0, nrow = 3, ncol = 2)
###############################
#   row  Category 1 #
###############################
for (i in 1:3){
compsR[i, 1] <- ZtZ[i, i]#Tube Location Component for Category 2#
compsR[i, 2] <- 1 - pchisq(compsR[i, 1], nc)#P-value of the Location Comp for Category 2#
}
compsRlast1 <- factor -sum (compsR[, 1])#Error of Row Components for Category 1#
if(nrow(Z) > 3) {
compsRlast2 <- 1 - pchisq(compsRlast1, (nr-3) * nc )#P-value for the Errors of Category 1#
Crow=rbind(compsR,c(compsRlast1,compsRlast2))
}
else {
compsRlast2 <- 0
Crow=rbind(compsR,c(compsRlast1,compsRlast2))
}
compsRtot1 <- factor
compsRtot2 <- 1 - pchisq(compsRtot1, nr * nc)
Crow=rbind(Crow,c(compsRtot1,compsRtot2))
###############################
#  column     Category 2 #
###############################
for (j in 1:3){
compsC[j, 1] <- tZZ[j, j]# Location Component for Category 1#
compsC[j, 2] <- 1 - pchisq(compsC[j, 1], nr)#P-Value for the Location Comp of Category 1#
}
compsClast1 <- factor -sum (compsC[, 1])#Error of Row Components for Category 1#
if(ncol(Z) > 3) {
compsClast2 <- 1 - pchisq(compsClast1, (nc-3) * nr )#P-value for the Errors of Category 1#
Ccol=rbind(compsC,c(compsClast1,compsClast2))
}
else {
compsClast2 <- 0
Ccol=rbind(compsC,c(compsClast1,compsClast2))
}
compsCtot1 <- factor
compsCtot2 <- 1 - pchisq(compsCtot1, nr * nc)
Ccol=rbind(Ccol,c(compsCtot1,compsCtot2))
list(compsR=Crow,compsC=Ccol)
#return(comps)

}
