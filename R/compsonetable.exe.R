compsonetable.exe <- 
function(Z){
ZtZ <- Z %*% t(Z)
factor <- sum(diag(ZtZ))
comps <- matrix(0, nrow=5, ncol=2)
###############################
#     Column Category         #
###############################
comps[2, 1] <- ZtZ[1, 1]#Location Component for Category 2#
comps[2, 2] <- 1 - pchisq(comps[2, 1], ncol(Z))#P-value of the Location Comp for Category 2#
comps[3, 1] <- ZtZ[2, 2]#Dispersion Component for Category 2#
comps[3, 2] <- 1 - pchisq(comps[3, 1], ncol(Z))#P-value of Dispersion Comp of Category 2#
comps[4, 1] <- factor - (comps[2, 1] + comps[3, 1])#Error of Components for Category 2#
if(nrow(Z) > 2) {
comps[4, 2] <- 1 - pchisq(comps[4, 1], (nrow(Z) - 2) * ncol(Z))#P-value of Errors for Category 2#
}
else {
comps[4, 2] <- 0
}
comps[5, 1] <- factor
comps[5, 2] <- 1 - pchisq(comps[5, 1], nrow(Z) * ncol(Z))
return(comps)
}
