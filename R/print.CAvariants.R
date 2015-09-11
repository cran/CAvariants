print.CAvariants <-
function(x,printdims=3,...) {
d <- min(printdims, x$cacorpo@S@r)
axnames <- character(length=d)
for (i in 1:d) { axnames[i] <- paste(" Axis",i) } 
cat("\n    RESULTS for",x$cacorpo@catype,  "Correspondence Analysis\n")
cat("\n    Data Table:\n")
print(x$cacorpo@DataMatrix)
cat("\n    Row Weights: Imass\n")
round(print(matrix(x$cacorpo@Imass,x$cacorpo@rows,x$cacorpo@rows,
dimnames=list(x$cacorpo@rowlabels,x$cacorpo@rowlabels)),digits=3))
cat("\n    Column Weights: Jmass\n")
round(print(matrix(x$cacorpo@Jmass,x$cacorpo@cols,x$cacorpo@cols,
dimnames=list(x$cacorpo@collabels,x$cacorpo@collabels)),digits=3))
cat("\n Total inertia ", round(x$cacorpo@inertiasum,digits=3), "\n\n")
cat("Inertias, percent inertias and cumulative percent inertias of the row space\n\n")
print(data.frame(x$cacorpo@inertias))
#----------------------------------------------------------------------------------------------
if ((x$cacorpo@catype=="DONSCA")|(x$cacorpo@catype=="SONSCA") ){
cat("Inertias, percent inertias and  cumulative percent inertias of the column space \n\n")
print(data.frame(x$cacorpo@inertias2))
cat("\n Diagonal element of the squared generalized correlation matrix:\n  Z'Z \n")
print(round(diag(t(x$cacorpo@Z)%*%x$cacorpo@Z),digits=6))
cat("\n Diagonal element of the squared generalized correlation matrix for a constant:\n  Z'Z for (n-1)*(I-1) \n")
print(round((sum(x$cacorpo@DataMatrix)-1)*(nrow(x$cacorpo@DataMatrix)-1)*diag(t(x$cacorpo@Z)%*%x$cacorpo@Z),digits=6))
cat("\n Diagonal element of the squared generalized correlation matrix:\n ZZ' \n")
print(round(diag(x$cacorpo@Z%*%t(x$cacorpo@Z)),digits=6))
cat("\n Diagonal element of the squared generalized correlation matrix for a constant:\n  ZZ' for (n-1)*(I-1) \n")
print(round((sum(x$cacorpo@DataMatrix)-1)*(nrow(x$cacorpo@DataMatrix)-1)*diag(x$cacorpo@Z%*%t(x$cacorpo@Z)),digits=6))
cat("\n Polynomial Components of Inertia \n")
print(x$cacorpo@comps)
}
#-----------------------------------------------------------------------------------------------
if ((x$cacorpo@catype=="DOCA")|(x$cacorpo@catype=="SOCA") ){
cat("Inertias, percent inertias and  cumulative percent inertias of the column space \n\n")
print(data.frame(x$cacorpo@inertias2))
cat("\n Diagonal element of the squared generalized correlation matrix:\n  Z'Z \n")
print(round(diag(t(x$cacorpo@Z)%*%x$cacorpo@Z),digits=6))
cat("\n Diagonal element of the squared generalized correlation matrix for a constant:\n  Z'Z for n \n")
print(round(diag(t(sum(x$cacorpo@DataMatrix)*x$cacorpo@Z)%*%x$cacorpo@Z),digits=6))
cat("\n Diagonal element of the squared generalized correlation matrix:\n ZZ' \n")
print(round(diag(x$cacorpo@Z%*%t(x$cacorpo@Z)),digits=6))
cat("\n Diagonal element of the squared generalized correlation matrix for a constant:\n  ZZ' for n \n")
print(round(diag(sum(x$cacorpo@DataMatrix)*x$cacorpo@Z%*%t(x$cacorpo@Z)),digits=6))
cat("\n Polynomial Components of Inertia \n")
print(x$cacorpo@comps)

}

#############################################################
if ((x$cacorpo@catype=="NSCA")||(x$cacorpo@catype=="DONSCA")||(x$cacorpo@catype=="SONSCA")){
cat("\n    Predictability Index for Variants of Non symmetrical Correspondence Analysis:\n")
cat("\nTau Index predicting from column \n\n")
print(x$cacorpo@S@tau)
Cstatistic<-(sum(x$cacorpo@DataMatrix)-1)*(nrow(x$cacorpo@DataMatrix)-1)*x$cacorpo@S@tau
pvalueC<-1 - pchisq(Cstatistic, (nrow(x$cacorpo@DataMatrix)-1)*(ncol(x$cacorpo@DataMatrix)-1))
cat("\n C-statistic", Cstatistic, "and p-value", pvalueC, "\n")
}
if ((x$cacorpo@catype=="DOCA")|(x$cacorpo@catype=="DONSCA")){
cat("\n Generalized correlation matrix of Bivariate Moment Decomposition\n")
print(round(x$cacorpo@Z,digits=3))
cat("\n Column polynomial axes \n")
printwithaxes(data.frame(x$cacorpo@S@Raxes[ ,1:d], row.names=x$cacorpo@collabels), axnames)
cat("\n Row polynomial axes \n")
printwithaxes(data.frame(x$cacorpo@S@Caxes[ ,1:d], row.names=x$cacorpo@rowlabels), axnames)
cat("\n Column standard polynomial coordinates \n")
printwithaxes(data.frame(x$cacorpo@Cstdcoord[, 1:d], row.names=x$cacorpo@collabels), axnames)
cat("\n Row standard polynomial coordinates \n")
printwithaxes(data.frame(x$cacorpo@Rstdcoord[, 1:d], row.names=x$cacorpo@rowlabels), axnames)
cat("\n Column principal polynomial coordinates \n")
printwithaxes(data.frame(x$cacorpo@Cprinccoord[, 1:d], row.names=x$cacorpo@collabels), axnames)
cat("\n Row principal polynomial coordinates \n")
printwithaxes(data.frame(x$cacorpo@Rprinccoord[, 1:d], row.names=x$cacorpo@rowlabels), axnames)
}

if ((x$cacorpo@catype=="SOCA")|(x$cacorpo@catype=="SONSCA")){
cat("\n Generalized correlation matrix of Hybrid Decomposition\n")
print(round(x$cacorpo@Z,digits=3))
cat("\n Row principal axes \n")
printwithaxes(data.frame(x$cacorpo@S@Caxes[ ,1:d], row.names=x$cacorpo@rowlabels), axnames)
cat("\n Column polynomial axes \n")
printwithaxes(data.frame(x$cacorpo@S@Raxes[ ,1:d], row.names=x$cacorpo@collabels), axnames)
cat("\n Column standard polynomial coordinates \n")
printwithaxes(data.frame(x$cacorpo@Cstdcoord[, 1:d], row.names=x$cacorpo@collabels), axnames)
cat("\n Row standard  coordinates \n")
printwithaxes(data.frame(x$cacorpo@Rstdcoord[, 1:d], row.names=x$cacorpo@rowlabels), axnames)
cat("\n Column principal  coordinates \n")
printwithaxes(data.frame(x$cacorpo@Cprinccoord[, 1:d], row.names=x$cacorpo@collabels), axnames)
cat("\n Row principal polynomial coordinates \n")
printwithaxes(data.frame(x$cacorpo@Rprinccoord[, 1:d], row.names=x$cacorpo@rowlabels), axnames)

}
else{
cat("\n Row principal axes \n")
printwithaxes(data.frame(x$cacorpo@S@Caxes[ ,1:d], row.names=x$cacorpo@rowlabels), axnames)
cat("\n Column principal axes \n")
printwithaxes(data.frame(x$cacorpo@S@Raxes[ ,1:d], row.names=x$cacorpo@collabels), axnames)
cat("\n Column standard coordinates \n")
printwithaxes(data.frame(x$cacorpo@Cstdcoord[, 1:d], row.names=x$cacorpo@collabels), axnames)
cat("\n Row standard coordinates \n")
printwithaxes(data.frame(x$cacorpo@Rstdcoord[, 1:d], row.names=x$cacorpo@rowlabels), axnames)
cat("\n Column principal  coordinates \n")
printwithaxes(data.frame(x$cacorpo@Cprinccoord[, 1:d], row.names=x$cacorpo@collabels), axnames)
cat("\n Row principal coordinates \n")
printwithaxes(data.frame(x$cacorpo@Rprinccoord[, 1:d], row.names=x$cacorpo@rowlabels), axnames)

}

cat("\n Column distances from the origin of the plot\n")
printwithaxes(data.frame((x$cacorpo@Cprinccoord^2)[, 1:d], row.names=x$cacorpo@collabels), axnames)
cat("\n Row distances from the origin of the plot \n")
printwithaxes(data.frame((x$cacorpo@Rprinccoord^2)[, 1:d], row.names=x$cacorpo@rowlabels), axnames)

cat("\n Inner product of coordinates (first two axes)   \n")
print(round(x$cacorpo@Trend,digits=3))

}
