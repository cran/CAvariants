printcacorporateplus <-
function(x) {
d <- min(x@printdims, x@S@r)
axnames <- character(length=d)
for (i in 1:d) { axnames[i] <- paste(" Axis",i) } 
cat("\n    RESULTS for",x@catype,  "Correspondence Analysis\n")
cat("\n    Data Table:\n")
print(x@DataMatrix)
cat("\n    Row Weights: Imass\n")
round(print(matrix(x@Imass,x@rows,x@rows,dimnames=list(x@rowlabels,x@rowlabels)),digits=3))
cat("\n    Column Weights: Jmass\n")
round(print(matrix(x@Jmass,x@cols,x@cols,dimnames=list(x@collabels,x@collabels)),digits=3))
cat("\n Total inertia ", round(x@inertiasum,digits=3), "\n\n")
cat("Inertias, percent inertias and cumulative percent inertias of the row space\n\n")
print(data.frame(x@inertias))
#----------------------------------------------------------------------------------------------
if ((x@catype=="DONSCA")|(x@catype=="SONSCA") ){
cat("Inertias, percent inertias and  cumulative percent inertias of the column space \n\n")
print(data.frame(x@inertias2))
cat("\n Diagonal element of the squared generalized correlation matrix:\n  Z'Z \n")
print(round(diag(t(x@Z)%*%x@Z),digits=6))
cat("\n Diagonal element of the squared generalized correlation matrix for a constant:\n  Z'Z for (n-1)*(I-1) \n")
print(round((sum(x@DataMatrix)-1)*(nrow(x@DataMatrix)-1)*diag(t(x@Z)%*%x@Z),digits=6))
cat("\n Diagonal element of the squared generalized correlation matrix:\n ZZ' \n")
print(round(diag(x@Z%*%t(x@Z)),digits=6))
cat("\n Diagonal element of the squared generalized correlation matrix for a constant:\n  ZZ' for (n-1)*(I-1) \n")
print(round((sum(x@DataMatrix)-1)*(nrow(x@DataMatrix)-1)*diag(x@Z%*%t(x@Z)),digits=6))
cat("\n Polynomial Components of Inertia \n")
print(x@comps)
}
#-----------------------------------------------------------------------------------------------
if ((x@catype=="DOCA")|(x@catype=="SOCA") ){
cat("Inertias, percent inertias and  cumulative percent inertias of the column space \n\n")
print(data.frame(x@inertias2))
cat("\n Diagonal element of the squared generalized correlation matrix:\n  Z'Z \n")
print(round(diag(t(x@Z)%*%x@Z),digits=6))
cat("\n Diagonal element of the squared generalized correlation matrix for a constant:\n  Z'Z for n \n")
print(round(diag(t(sum(x@DataMatrix)*x@Z)%*%x@Z),digits=6))
cat("\n Diagonal element of the squared generalized correlation matrix:\n ZZ' \n")
print(round(diag(x@Z%*%t(x@Z)),digits=6))
cat("\n Diagonal element of the squared generalized correlation matrix for a constant:\n  ZZ' for n \n")
print(round(diag(sum(x@DataMatrix)*x@Z%*%t(x@Z)),digits=6))
cat("\n Polynomial Components of Inertia \n")
print(x@comps)

}

#############################################################
if ((x@catype=="NSCA")||(x@catype=="DONSCA")||(x@catype=="SONSCA")){
cat("\n    Predictability Index for Variants of Non symmetrical Correspondence Analysis:\n")
cat("\nTau Index predicting from column \n\n")
print(x@S@tau)
Cstatistic<-(sum(x@DataMatrix)-1)*(nrow(x@DataMatrix)-1)*x@S@tau
pvalueC<-1 - pchisq(Cstatistic, (nrow(x@DataMatrix)-1)*(ncol(x@DataMatrix)-1))
cat("\n C-statistic", Cstatistic, "and p-value", pvalueC, "\n")
}
if ((x@catype=="DOCA")|(x@catype=="DONSCA")){
cat("\n Generalized correlation matrix of Bivariate Moment Decomposition\n")
print(round(x@Z,digits=3))
cat("\n Column polynomial axes \n")
printwithaxes(data.frame(x@S@Raxes[ ,1:d], row.names=x@collabels), axnames)
cat("\n Row polynomial axes \n")
printwithaxes(data.frame(x@S@Caxes[ ,1:d], row.names=x@rowlabels), axnames)
cat("\n Column standard polynomial coordinates \n")
printwithaxes(data.frame(x@Cstdcoord[, 1:d], row.names=x@collabels), axnames)
cat("\n Row standard polynomial coordinates \n")
printwithaxes(data.frame(x@Rstdcoord[, 1:d], row.names=x@rowlabels), axnames)
cat("\n Column principal polynomial coordinates \n")
printwithaxes(data.frame(x@Cprinccoord[, 1:d], row.names=x@collabels), axnames)
cat("\n Row principal polynomial coordinates \n")
printwithaxes(data.frame(x@Rprinccoord[, 1:d], row.names=x@rowlabels), axnames)

}

if ((x@catype=="SOCA")|(x@catype=="SONSCA")){
cat("\n Generalized correlation matrix of Hybrid Decomposition\n")
print(round(x@Z,digits=3))
cat("\n Row principal axes \n")
printwithaxes(data.frame(x@S@Caxes[ ,1:d], row.names=x@rowlabels), axnames)
cat("\n Column polynomial axes \n")
printwithaxes(data.frame(x@S@Raxes[ ,1:d], row.names=x@collabels), axnames)
cat("\n Column standard polynomial coordinates \n")
printwithaxes(data.frame(x@Cstdcoord[, 1:d], row.names=x@collabels), axnames)
cat("\n Row standard  coordinates \n")
printwithaxes(data.frame(x@Rstdcoord[, 1:d], row.names=x@rowlabels), axnames)
cat("\n Column principal  coordinates \n")
printwithaxes(data.frame(x@Cprinccoord[, 1:d], row.names=x@collabels), axnames)
cat("\n Row principal polynomial coordinates \n")
printwithaxes(data.frame(x@Rprinccoord[, 1:d], row.names=x@rowlabels), axnames)

}
else{
cat("\n Row principal axes \n")
printwithaxes(data.frame(x@S@Caxes[ ,1:d], row.names=x@rowlabels), axnames)
cat("\n Column principal axes \n")
printwithaxes(data.frame(x@S@Raxes[ ,1:d], row.names=x@collabels), axnames)
cat("\n Column standard coordinates \n")
printwithaxes(data.frame(x@Cstdcoord[, 1:d], row.names=x@collabels), axnames)
cat("\n Row standard coordinates \n")
printwithaxes(data.frame(x@Rstdcoord[, 1:d], row.names=x@rowlabels), axnames)
cat("\n Column principal  coordinates \n")
printwithaxes(data.frame(x@Cprinccoord[, 1:d], row.names=x@collabels), axnames)
cat("\n Row principal coordinates \n")
printwithaxes(data.frame(x@Rprinccoord[, 1:d], row.names=x@rowlabels), axnames)

}


cat("\n Inner product of coordinates (first two axes)   \n")
print(round(x@Trend,digits=3))

}
