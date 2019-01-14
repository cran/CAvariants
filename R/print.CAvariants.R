print.CAvariants <-
function(x,printdims=2,ellprint=TRUE, alpha=0.05,digits=3,...) {
#d <- min(printdims, x$r)
d<-printdims
if (d>x$r) { d<-x$r
cat("The maximum dimension number cannot be greater than the rank of the matrix\n")}
axnames <- character(length=d)
for (i in 1:d) { axnames[i] <- paste(" Axis",i) } 
cat("\n    RESULTS for",x$catype,  "Correspondence Analysis\n")
cat("\n    Data Matrix:\n")
print(x$Xtable)
cat("\n    Row weights: Imass\n")
print(round(matrix(x$Imass,x$rows,x$rows,
dimnames=list(x$rowlabels,x$rowlabels)),digits=digits))
cat("\n    Column Weights: Jmass\n")
print(round(matrix(x$Jmass,x$cols,x$cols,
dimnames=list(x$collabels,x$collabels)),digits=digits))
#---------------------------------------------------------------------------
if ((x$catype=="CA")|(x$catype=="NSCA") ){
cat("\n Total inertia ", round(x$inertiasum,digits=digits), "\n\n")
cat("The inertia values, their percentage contribution to the total inertia and   
the cumulative percent inertias  \n")
print(round(data.frame(x$inertias),digits=digits))
}
#----------------------------------------------------------------------------------------------
if ((x$catype=="DONSCA")|(x$catype=="SONSCA") ){
cat("\n Total inertia ", round(x$inertiasum,digits=digits), "\n\n")
cat("Inertias, percent inertias and cumulative percent inertias of the row  space\n\n")
print(round(data.frame(x$inertias),digits=digits))
cat("Inertias, percent inertias and  cumulative percent inertias of the column space \n\n")
print(round(data.frame(x$inertias2),digits=digits))
cat("\n Absolute Contributions of Rows (per 100):\n")
ctrRow<-x$Imass%*%x$Rprinccoord^2%*%diag(1/x$inertias[,1])*100
print(round(ctrRow,digits=digits))
cat("\n Absolute Contributions of Columns (per 100):\n")
ctrCol<-x$Jmass%*%x$Cprinccoord^2%*%diag(1/x$inertias[,1])*100
print(round(ctrCol,digits=digits))
cat("\n Relative Contributions of Rows (per 100):\n")
ctrRowrel<-x$Rprinccoord^2/apply(x$Rprinccoord^2,1,sum)*100
printwithaxes(data.frame((ctrRowrel)[, 1:d], row.names=x$rowlabels), axnames,digits=digits)
cat("\n Relative Contributions of Columns (per 100):\n")
ctrColrel<-x$Cprinccoord^2/apply(x$Cprinccoord^2,1,sum)*100
printwithaxes(data.frame((ctrColrel)[, 1:d], row.names=x$collabels), axnames,digits=digits)
}
#-----------------------------------------------------------------------------------------------
if ((x$catype=="DOCA")|(x$catype=="SOCA") ){
cat("\n Total inertia ", round(x$inertiasum,digits=digits), "\n\n")
cat("Inertias, percent inertias and cumulative percent inertias of the row space\n\n")
print(round(data.frame(x$inertias),digits=digits))
ctrRow<-x$Imass%*%x$Rprinccoord^2%*%diag(1/x$inertias[,1])*100
cat("Inertias, percent inertias and  cumulative percent inertias of the column space \n\n")
print(round(data.frame(x$inertias2),digits=digits))
ctrCol<-x$Jmass%*%x$Cprinccoord^2%*%diag(1/x$inertias[,1])*100
cat("\n Absolute Contributions of Columns (per 100):\n")
print(round(ctrCol,digits=digits))
cat("\n Absolute Contributions of Rows (per 100):\n")
print(round(ctrRow,digits=digits))
cat("\n Relative Contributions of Rows (per 100):\n")
ctrRowrel<-x$Rprinccoord^2/apply(x$Rprinccoord^2,1,sum)*100
printwithaxes(data.frame((ctrRowrel)[, 1:d], row.names=x$rowlabels), axnames,digits=digits)
cat("\n Relative Contributions of Columns (per 100):\n")
ctrColrel<-x$Cprinccoord^2/apply(x$Cprinccoord^2,1,sum)*100
printwithaxes(data.frame((ctrColrel)[, 1:d], row.names=x$collabels), axnames,digits=digits)
}
#############################################################
if ((x$catype=="NSCA")||(x$catype=="DONSCA")||(x$catype=="SONSCA")){
cat("\n    Predictability Index for Variants of Non symmetrical Correspondence Analysis:\n")
cat("\n Numerator of Tau Index predicting the row categories from the column categories\n\n")
print(round(x$inertiasum,digits=digits))
cat("\n Tau Index predicting the row categories from the column categories\n\n")
print(round(x$inertiasum/x$tauden,digits=digits))
Cstatistic<-(sum(x$Xtable)-1)*(nrow(x$Xtable)-1)*x$tau
pvalueC<-1 - pchisq(Cstatistic, (nrow(x$Xtable)-1)*(ncol(x$Xtable)-1))
cat("\n C-statistic", Cstatistic, "and p-value", pvalueC, "\n")
}
if ((x$catype=="DOCA")|(x$catype=="DONSCA")){
cat("\n Polynomial Components of Inertia \n
** Row Components ** \n")
print(round(x$comps$compsR,digits=digits))
cat("\n** Column Components ** \n")
print(round(x$comps$compsC,digits=digits))
cat("\n Generalized correlation matrix of Bivariate Moment Decomposition\n")
print(round(x$Z,digits=digits))
#cat("\n Column polynomial axes \n")
#printwithaxes(data.frame(x$Raxes[ ,1:d], row.names=x$collabels), axnames)
#cat("\n Row polynomial axes \n")
#printwithaxes(data.frame(x$Caxes[ ,1:d], row.names=x$rowlabels), axnames)
cat("\n Column standard polynomial coordinates = column polynomial axes \n")
printwithaxes(data.frame(x$Cstdcoord[, 1:d], row.names=x$collabels), axnames,digits=digits)
cat("\n Row standard polynomial coordinates = row polynomial axes  \n")
printwithaxes(data.frame(x$Rstdcoord[, 1:d], row.names=x$rowlabels), axnames,digits=digits)
cat("\n Column principal polynomial coordinates \n")
printwithaxes(data.frame(x$Cprinccoord[, 1:d], row.names=x$collabels), axnames,digits=digits)
cat("\n Row principal polynomial coordinates \n")
printwithaxes(data.frame(x$Rprinccoord[, 1:d], row.names=x$rowlabels), axnames,digits=digits)
}

if ((x$catype=="SOCA")|(x$catype=="SONSCA")){
cat("\n Polynomial Components of Inertia \n
** Column Components ** \n")
print(round(x$comps$comps,digits=digits))
cat("\n Generalized correlation matrix of Hybrid Decomposition\n")
print(round(x$Z,digits=digits))
cat("\n Column standard polynomial coordinates  \n")
printwithaxes(data.frame(x$Cstdcoord[, 1:d], row.names=x$collabels), axnames,digits=digits)
#cat("\n Row standard  coordinates  = row principal axes\n")
#printwithaxes(data.frame(x$Rstdcoord[, 1:d], row.names=x$rowlabels), axnames,digits=digits)
#cat("\n Column principal  coordinates \n")
#printwithaxes(data.frame(x$Cprinccoord[, 1:d], row.names=x$collabels), axnames,digits=digits)
cat("\n Row principal polynomial coordinates \n")
printwithaxes(data.frame(x$Rprinccoord[, 1:d], row.names=x$rowlabels), axnames,digits=digits)
}
if ((x$catype=="CA")|(x$catype=="NSCA")){
cat("\n Column standard coordinates \n")
printwithaxes(data.frame(x$Cstdcoord[, 1:d], row.names=x$collabels), axnames,digits=digits)
cat("\n Row standard coordinates \n")
printwithaxes(data.frame(x$Rstdcoord[, 1:d], row.names=x$rowlabels), axnames,digits=digits)
cat("\n Column principal  coordinates \n")
printwithaxes(data.frame(x$Cprinccoord[, 1:d], row.names=x$collabels), axnames,digits=digits)
cat("\n Row principal coordinates \n")
printwithaxes(data.frame(x$Rprinccoord[, 1:d], row.names=x$rowlabels), axnames,digits=digits)
cat("\n Absolute Contributions of Rows (per 100):\n")
ctrRow<-x$Imass%*%x$Rprinccoord^2%*%diag(1/x$inertias[,1])*100
print(round(ctrRow,digits=digits))
cat("\n Absolute Contributions of Columns (per 100):\n")
ctrCol<-x$Jmass%*%x$Cprinccoord^2%*%diag(1/x$inertias[,1])*100
print(round(ctrCol,digits=digits))
cat("\n Relative Contributions of Rows (per 100):\n")
ctrRowrel<-x$Rprinccoord^2/apply(x$Rprinccoord^2,1,sum)*100
printwithaxes(data.frame((ctrRowrel^2)[, 1:d], row.names=x$rowlabels), axnames,digits=digits)
cat("\n Relative Contributions of Columns (per 100):\n")
ctrColrel<-x$Cprinccoord^2/apply(x$Cprinccoord^2,1,sum)*100
printwithaxes(data.frame((ctrColrel)[, 1:d], row.names=x$collabels), axnames,digits=digits)

}
cat("\n Column distances from the origin of the plot\n")
printwithaxes(data.frame((x$Cprinccoord^2)[, 1:d], row.names=x$collabels), axnames,digits=digits)
cat("\n Row distances from the origin of the plot \n")
printwithaxes(data.frame((x$Rprinccoord^2)[, 1:d], row.names=x$rowlabels), axnames,digits=digits)
cat("\n Inner product of coordinates (first two axes when 'firstaxis=1' and 'lastaxis=2')   \n")
print(round(x$Trend,digits=digits))
#browser()
if (ellprint==TRUE){
#dimnames(x$resellprint$row.summ)[[1]]<-dimnames(x$Xtable)[[1]]
#dimnames(x$resellprint$col.summ)[[1]]<-dimnames(x$Xtable)[[2]]
cat("\n    Eccentricity of ellipses\n")
print(round(x$resellprint$eccentricity,digits=digits))
cat("\n    Ellipse axes, Area, p-values of rows\n")
print(round(x$resellprint$row.summ,digits=digits))
cat("\n    Ellipse axes, Area, p-values of columns\n")
print(round(x$resellprint$col.summ,digits=digits))
}#end ell

}
