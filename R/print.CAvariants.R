print.CAvariants <-
function(x,printdims=2,ellprint=TRUE,Mell=min(nrow(x$Xtable),ncol(x$Xtable))-1,alpha=0.05,digits=3,...) {
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
#cat("\n Diagonal element of the squared generalized correlation matrix:\n  Z'Z \n")
#print(round(diag(t(x$Z)%*%x$Z),digits=digits))
#cat("\n Diagonal element of the squared generalized correlation matrix for a constant:\n  Z'Z for (n-1)*(I-1) \n")
#print(round((sum(x$Xtable)-1)*(nrow(x$Xtable)-1)*diag(t(x$Z)%*%x$Z),digits=digits))
#cat("\n Diagonal element of the squared generalized correlation matrix:\n ZZ' \n")
#print(round(diag(x$Z%*%t(x$Z)),digits=digits))
#cat("\n Diagonal element of the squared generalized correlation matrix for a constant:\n  ZZ' for (n-1)*(I-1) \n")
#print(round((sum(x$Xtable)-1)*(nrow(x$Xtable)-1)*diag(x$Z%*%t(x$Z)),digits=digits))
#cat("\n Polynomial Components of Inertia \n")
#print(x$comps)
}
#-----------------------------------------------------------------------------------------------
if ((x$catype=="DOCA")|(x$catype=="SOCA") ){
cat("\n Total inertia ", round(x$inertiasum,digits=digits), "\n\n")
cat("Inertias, percent inertias and cumulative percent inertias of the row space\n\n")
print(round(data.frame(x$inertias),digits=digits))
cat("Inertias, percent inertias and  cumulative percent inertias of the column space \n\n")
print(round(data.frame(x$inertias2),digits=digits))
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
cat("\n Column standard polynomial coordinates = column polynomial axes \n")
printwithaxes(data.frame(x$Cstdcoord[, 1:d], row.names=x$collabels), axnames,digits=digits)
#cat("\n Row standard  coordinates  = row principal axes\n")
#printwithaxes(data.frame(x$Rstdcoord[, 1:d], row.names=x$rowlabels), axnames,digits=digits)
#cat("\n Column principal  coordinates \n")
#printwithaxes(data.frame(x$Cprinccoord[, 1:d], row.names=x$collabels), axnames,digits=digits)
cat("\n Row principal polynomial coordinates \n")
printwithaxes(data.frame(x$Rprinccoord[, 1:d], row.names=x$rowlabels), axnames,digits=digits)
}
if ((x$catype=="CA")|(x$catype=="NSCA")){
cat("\n Column standard coordinates = column principal axes\n")
printwithaxes(data.frame(x$Cstdcoord[, 1:d], row.names=x$collabels), axnames,digits=digits)
cat("\n Row standard coordinates = row principal axes \n")
printwithaxes(data.frame(x$Rstdcoord[, 1:d], row.names=x$rowlabels), axnames,digits=digits)
cat("\n Column principal  coordinates \n")
printwithaxes(data.frame(x$Cprinccoord[, 1:d], row.names=x$collabels), axnames,digits=digits)
cat("\n Row principal coordinates \n")
printwithaxes(data.frame(x$Rprinccoord[, 1:d], row.names=x$rowlabels), axnames,digits=digits)
}
cat("\n Column distances from the origin of the plot\n")
printwithaxes(data.frame((x$Cprinccoord^2)[, 1:d], row.names=x$collabels), axnames,digits=digits)
cat("\n Row distances from the origin of the plot \n")
printwithaxes(data.frame((x$Rprinccoord^2)[, 1:d], row.names=x$rowlabels), axnames,digits=digits)
cat("\n Inner product of coordinates (first two axes when 'firstaxis=1' and 'lastaxis=2')   \n")
print(round(x$Trend,digits=digits))

if (ellprint==TRUE){

cord1<-x$Cprinccoord
cord2<-x$Rprinccoord
if ((x$catype=="DOCA")|(x$catype=="SOCA")|(x$catype=="SONSCA")|(x$catype=="DONSCA") ){
cordr<-cord2
cordc<-cord1
cord1<-cordr
cord2<-cordc
}
risell<-switch(x$catype, "CA"=caellipse(Xtable=x$Xtable,a1=1,a2=2,alpha=alpha,M=Mell,prop=1,
Imass=x$Imass,Jmass=x$Jmass,a=x$Rstdcoord,b=x$Cstdcoord,g=cord1,fr=cord2,dmu=diag(x$inertias[,1]),inertiapc=x$inertias[,2],
plottype="biplot",biptype="row",pos=1,arrow=T,length=0,graphy=F), 

"SOCA"=caellipse(Xtable=x$Xtable,a1=1,a2=2,alpha=alpha,M=Mell,prop=1,
Imass=x$Imass,Jmass=x$Jmass,a=x$Rstdcoord,b=x$Cstdcoord,g=cord2,fr=cord1,dmu=diag(x$inertias2[,1]),inertiapc=x$inertias2[,2],
plottype="biplot",biptype="row",pos=1,arrow=T,length=0,graphy=F), 

"DOCA"=caellipse(Xtable=x$Xtable,a1=1,a2=2,alpha=alpha,M=Mell,prop=1,
Imass=x$Imass,Jmass=x$Jmass,a=x$Rstdcoord,
b=x$Cstdcoord,g=cord2,fr=cord1,dmu=diag(x$inertias2[,1]),inertiapc=x$inertias2[,2],
plottype="biplot",biptype="row",pos=1,arrow=T,length=0,graphy=F), 

"NSCA"=nscaellipse(Xtable=x$Xtable,a1=1,a2=2,alpha=alpha,M=Mell,prop=1,
Imass=x$Imass,Jmass=x$Jmass,a=x$Rstdcoord,
b=x$Cstdcoord,g=cord1,fr=cord2,dmu=diag(x$inertias[,1]), tauden=x$tauden,
inertiapc=x$inertias[,2],
plottype="biplot",biptype="row",pos=1,arrow=T,length=0,graphy=F),

"SONSCA"=nscaellipse(Xtable=x$Xtable,a1=1,a2=2,alpha=alpha,M= Mell,prop=1,
Imass=x$Imass,Jmass=x$Jmass,a=x$Rstdcoord,
b=x$Cstdcoord,g=cord2,fr=cord1,dmu=diag(x$inertias2[,1]),tauden=x$tauden,inertiapc=x$inertias2[,2],
plottype="biplot",biptype="row",pos=1,arrow=T,length=0,graphy=F), 

"DONSCA"=nscaellipse(Xtable=x$Xtable,a1=1,a2=2,alpha=alpha,M=Mell,prop=1,
Imass=x$Imass,Jmass=x$Jmass,a=x$Rstdcoord,
b=x$Cstdcoord,g=cord2,fr=cord1,dmu=diag(x$inertias2[,1]),tauden=x$tauden,inertiapc=x$inertias2[,2],
plottype="biplot",biptype="row",pos=1,arrow=F,length=0,graphy=F))
dimnames(risell$row.summ)[[1]]<-dimnames(x$Xtable)[[1]]
dimnames(risell$col.summ)[[1]]<-dimnames(x$Xtable)[[2]]
cat("\n    Eccentricity of ellipses\n")
print(round(risell$eccentricity,digits=digits))
cat("\n    Ellipse axes, Area, p-values of rows\n")
print(round(risell$row.summ,digits=digits))
#print(round(risell$row.summ))
cat("\n    Ellipse axes, Area, p-values of columns\n")
print(round(risell$col.summ,digits=digits))
}#end ell

}
