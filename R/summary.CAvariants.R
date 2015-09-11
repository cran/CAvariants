summary.CAvariants <-
function(object,...) {
cat("\n    RESULTS for",object$cacorpo@catype,  "Correspondence Analysis\n")
cat("\n Total inertia ", round(object$cacorpo@inertiasum,digits=3), "\n\n")
cat("Inertias, percent inertias and cumulative percent inertias of the row space\n\n")
print(data.frame(object$cacorpo@inertias))
#----------------------------------------------------------------------------------------------
if ((object$cacorpo@catype=="DONSCA")|(object$cacorpo@catype=="SONSCA") ){
cat("Inertias, percent inertias and  cumulative percent inertias of the column space \n\n")
print(data.frame(object$cacorpo@inertias2))
cat("\n Polynomial Components of Inertia \n")
print(object$cacorpo@comps)
}
#-----------------------------------------------------------------------------------------------
if ((object$cacorpo@catype=="DOCA")|(object$cacorpo@catype=="SOCA") ){
cat("Inertias, percent inertias and  cumulative percent inertias of the column space \n\n")
print(data.frame(object$cacorpo@inertias2))
cat("\n Polynomial Components of Inertia \n")
print(object$cacorpo@comps)
}

#############################################################
if ((object$cacorpo@catype=="NSCA")||(object$cacorpo@catype=="DONSCA")||(object$cacorpo@catype=="SONSCA")){
cat("\n    Predictability Index for Variants of Non symmetrical Correspondence Analysis:\n")
cat("\nTau Index predicting from column \n\n")
print(object$cacorpo@S@tau)
Cstatistic<-(sum(object$cacorpo@DataMatrix)-1)*(nrow(object$cacorpo@DataMatrix)-1)*object$cacorpo@S@tau
pvalueC<-1 - pchisq(Cstatistic, (nrow(object$cacorpo@DataMatrix)-1)*(ncol(object$cacorpo@DataMatrix)-1))
cat("\n C-statistic", Cstatistic, "and p-value", pvalueC, "\n")
}
if ((object$cacorpo@catype=="DOCA")|(object$cacorpo@catype=="DONSCA")){
print(data.frame(object$cacorpo@S@Caxes, row.names=object$cacorpo@rowlabels), digits=3)
cat("\n Column standard polynomial coordinates \n")
print(data.frame(object$cacorpo@Cstdcoord, row.names=object$cacorpo@collabels), digits=3)
cat("\n Row standard polynomial coordinates \n")
print(data.frame(object$cacorpo@Rstdcoord, row.names=object$cacorpo@rowlabels), digits=3)
cat("\n Column principal polynomial coordinates \n")
print(data.frame(object$cacorpo@Cprinccoord, row.names=object$cacorpo@collabels), digits=3)
cat("\n Row principal polynomial coordinates \n")
print(data.frame(object$cacorpo@Rprinccoord, row.names=object$cacorpo@rowlabels), digits=3)
}

if ((object$cacorpo@catype=="SOCA")|(object$cacorpo@catype=="SONSCA")){
cat("\n Column standard polynomial coordinates \n")
print(data.frame(object$cacorpo@Cstdcoord, row.names=object$cacorpo@collabels), digits=3)
cat("\n Row standard  coordinates \n")
print(data.frame(object$cacorpo@Rstdcoord, row.names=object$cacorpo@rowlabels), digits=3)
cat("\n Column principal  coordinates \n")
print(data.frame(object$cacorpo@Cprinccoord, row.names=object$cacorpo@collabels), digits=3)
cat("\n Row principal polynomial coordinates \n")
print(data.frame(object$cacorpo@Rprinccoord, row.names=object$cacorpo@rowlabels), digits=3)

}
else{
cat("\n Row principal axes \n")
cat("\n Column standard coordinates \n")
print(data.frame(object$cacorpo@Cstdcoord, row.names=object$cacorpo@collabels), digits=3)
cat("\n Row standard coordinates \n")
print(data.frame(object$cacorpo@Rstdcoord, row.names=object$cacorpo@rowlabels), digits=3)
cat("\n Column principal  coordinates \n")
print(data.frame(object$cacorpo@Cprinccoord, row.names=object$cacorpo@collabels), digits=3)
cat("\n Row principal coordinates \n")
print(data.frame(object$cacorpo@Rprinccoord, row.names=object$cacorpo@rowlabels), digits=3)

}


cat("\n Inner product of coordinates (first two axes)   \n")
print(round(object$cacorpo@Trend,digits=3))

}
