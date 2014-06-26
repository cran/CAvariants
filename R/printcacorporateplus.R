printcacorporateplus <-
function(x) {
d <- min(x@printdims, x@S@r)
axnames <- character(length=d)
for (i in 1:d) { axnames[i] <- paste(" Axis",i) } 
cat("\n    RESULTS for",x@catype,  "Correspondence Analysis\n")
cat("\n    Data Table:\n")
print(x@DataMatrix)
cat("\n    Row Weights: Imass:\n")
round(print(matrix(x@Imass,x@rows,x@rows,dimnames=list(x@rowlabels,x@rowlabels)),digits=3))
cat("\n    Column Weights: Jmass:\n")
round(print(matrix(x@Jmass,x@cols,x@cols,dimnames=list(x@collabels,x@collabels)),digits=3))

cat("\n Total inertia ", round(x@inertiasum,digits=3), "\n\n")
cat("Inertias, percent inertias and cumulative percent inertias of the row space\n\n")
print(data.frame(x@inertias))
if ((x@catype=="DOCA")|(x@catype=="DONSCA")){
cat("Inertias, percent inertias and cumulative percent inertias of the column space \n\n")
print(data.frame(x@inertias2))}

cat("\nRow standard coordinates\n\n")
printwithaxes(data.frame(x@Rstdcoord[ ,1:d], row.names=x@rowlabels), axnames)
cat("\nRow principal coordinates\n\n")
printwithaxes(data.frame(x@Rprinccoord[ ,1:d], row.names=x@rowlabels), axnames)
cat("\nColumn standard coordinates\n\n")
printwithaxes(data.frame(x@Cstdcoord[ ,1:d], row.names=x@collabels), axnames)
cat("\nColumn principal coordinates\n\n")
printwithaxes(data.frame(x@Cprinccoord[ ,1:d], row.names=x@collabels), axnames)

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
cat("\n Column polynomial axes \n")
printwithaxes(data.frame(x@S@Raxes[ ,1:d], row.names=x@collabels), axnames)
cat("\n Row polynomial axes \n")
printwithaxes(data.frame(x@S@Caxes[ ,1:d], row.names=x@rowlabels), axnames)
}

if ((x@catype=="SOCA")|(x@catype=="SONSCA")){
cat("\n Row principal axes \n")
printwithaxes(data.frame(x@S@Caxes[ ,1:d], row.names=x@rowlabels), axnames)
cat("\n Column polynomial axes \n")
printwithaxes(data.frame(x@S@Raxes[ ,1:d], row.names=x@collabels), axnames)
}
#else{
#cat("\n Row principal axes \n")
#printwithaxes(data.frame(x@S@Caxes[ ,1:d], row.names=x@rowlabels), axnames)
#cat("\n Column principal axes \n")
#printwithaxes(data.frame(x@S@Raxes[ ,1:d], row.names=x@collabels), axnames)
#}
####################################################################
if((x@catype=="DOCA")||(x@catype=="SOCA")||(x@catype=="DONSCA")||(x@catype=="SONSCA")){
cat("\n Generalized correlation matrix \n")
print(round(x@Z,digits=3))
cat("\n Polynomial Components of Inertia \n")
print(x@comps)
}
else
{cat("\n\n")}

#if ((x@catype=="DOCA")|(x@catype=="SOCA"))
#{
#Trend<-(x@Rprinccoord[,1:2]%*%t(x@Jmass%*%x@Cstdcoord[,1:2]))
#}
#browser()
#if ((x@catype=="CA")|(x@catype=="NSCA")|(x@catype=="SONSCA")|(x@catype=="DONSCA")){

#Trend<-t(x@Cprinccoord[,1:2]%*%t(x@Imass%*%x@Rstdcoord[,1:2]))
#}
#dimnames(Trend)<-list(x@rowlabels,x@collabels)
cat("\n Inner product of coordinates (first two axes)   \n")
print(round(x@Trend,digits=3))

}
