compstable.exe <-
function(Z){
tZZ <- t(Z) %*% Z
ZtZ <- Z %*% t(Z)
factor <- sum(diag(tZZ))
comps <- matrix(0, nrow = 9, ncol = 2)
###############################
#     Category 1 #
###############################
comps[2, 1] <- tZZ[1, 1]# Location Component for Category 1#
comps[2, 2] <- 1 - pchisq(comps[2, 1], nrow(Z))#P-Value for the Location Comp of Category 1#
comps[3, 1] <- tZZ[2, 2]#Dispersion Component for Category 1#
comps[3, 2] <- 1 - pchisq(comps[3, 1], nrow(Z))#P-value for the Dispersion Comp of Category 2#
comps[4, 1] <- factor - (comps[2, 1] + comps[3, 1])#Error of Row Components for Category 1#
if(ncol(Z) > 2) {
comps[4, 2] <- 1 - pchisq(comps[4, 1], nrow(Z) * (ncol(Z) - 2))#P-value for the Errors of Category 1#
}
else {
comps[4, 2] <- 0
}
###############################
#     Category 2 #
###############################
comps[6, 1] <- ZtZ[1, 1]#Tube Location Component for Category 2#
comps[6, 2] <- 1 - pchisq(comps[6, 1], ncol(Z))#P-value of the Location Comp for Category 2#
comps[7, 1] <- ZtZ[2, 2]#Dispersion Component for Category 2#
comps[7, 2] <- 1 - pchisq(comps[7, 1], ncol(Z))#P-value of Dispersion Comp of Category 2#
comps[8, 1] <- factor - (comps[6, 1] + comps[7, 1])#Error of Components for Category 2#
if(nrow(Z) > 2) {
comps[8, 2] <- 1 - pchisq(comps[8, 1], (nrow(Z) - 2) * ncol(Z))#P-value of Errors for Category 2#
}
else {
comps[8, 2] <- 0
}
comps[9, 1] <- factor
comps[9, 2] <- 1 - pchisq(comps[9, 1], nrow(Z) * ncol(Z))
return(comps)
}
