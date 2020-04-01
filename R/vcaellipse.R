#library(ggplot2)
#library(ggforce)
#library(plotly)
#library(CAvariants)
#vcaellipse<-
#function(t.inertia=t.inertia,inertias=inertias[,1],inertiapc=#x$inertias[,2],cord1=x$Rprinccoord,cord2=x$Cprinccoord,a=x$Rstdcoord,b=x$Cstdcoord,firstaxis=1,lastaxis=2,n=x$n,M=min(nrow(x$Xtable),ncol(x$Xtable))-1,Imass=x#$Imass,Jmass=x$Jmass){
#-----------------
vcaellipse<-function(t.inertia,inertias,inertiapc,cord1,cord2,a,b,firstaxis=1,lastaxis=2,n,M=2,Imass,Jmass){
#globalVariables(categ)
Inames<-thinglabels<-dimnames(cord1)[[1]]
Jnames<-varlabels<-dimnames(cord2)[[1]]
nthings<-I<-dim(cord1)[1]
nvars<-J<-dim(cord2)[1]
dmu<-sqrt(inertias)
alpha=0.05
#-----------------------------axis ellipses
chisq.val <- qchisq(1 - alpha, df = (I - 1) * (J - 1))
hlax1.row <- vector(mode = "numeric", length = I)
hlax2.row <- vector(mode = "numeric", length = I)
hlax1.col <- vector(mode = "numeric", length = J)
hlax2.col <- vector(mode = "numeric", length = J)
#browser()
if (M > 2) {
  for (i in 1:I) {
    hlax1.row[i] <- dmu[1] * sqrt(abs((chisq.val/(t.inertia)) * 
                                          (1/Imass[i,i] - sum(a[i, 3:M]^2))))
    hlax2.row[i] <- dmu[2] * sqrt(abs((chisq.val/(t.inertia)) * 
                                          (1/Imass[i,i] - sum(a[i, 3:M]^2))))
  }
  for (j in 1:J) {
    hlax1.col[j] <- dmu[1] * sqrt(abs((chisq.val/(t.inertia)) * 
                                          (1/Jmass[j,j] - sum(b[j, 3:M]^2))))
    hlax2.col[j] <- dmu[2] * sqrt(abs((chisq.val/(t.inertia)) * 
                                          (1/Jmass[j,j] - sum(b[j, 3:M]^2))))
      }
            }
if (M == 2) {
    for (i in 1:I) {
    hlax1.row[i] <- dmu[1] * sqrt(abs((chisq.val/(t.inertia)) * 
                                          (1/(Imass)[i,i])))
    hlax2.row[i] <- dmu[2] * sqrt(abs((chisq.val/(t.inertia)) * 
                                          (1/(Imass)[i,i])))
  }
  for (j in 1:J) {
    hlax1.col[j] <- dmu[1] * sqrt(abs((chisq.val/(t.inertia)) * 
                                          (1/(Jmass)[j,j])))
    hlax2.col[j] <- dmu[2] * sqrt(abs((chisq.val/(t.inertia)) * 
                                          (1/(Jmass)[j,j])))
  }
}
#---------------------------------------------------
#faxis <- data.frame(coord=c(hlax1.row,hlax2.row), labels=thinglabels, categ=rep("rows", nthings)) # build a dataframe to be used as input for plotting via ggplot2
#gaxis <- data.frame(coord=c(hlax1.col,hlax2.col), labels=varlabels, categ=rep("cols", nvars)) # build a dataframe to be used as input for plotting via ggplot2
faxis <- cbind(hlax1.row,hlax2.row)
gaxis <- cbind(hlax1.col,hlax2.col)
categ<-NULL
frows <- data.frame(coord=cord1, labels=thinglabels, categ=rep("rows", nthings)) # build a dataframe to be used as input for plotting via ggplot2
gcols <- data.frame(coord=cord2, labels=varlabels, categ=rep("cols", nvars)) # build a dataframe to be used as input for plotting via ggplot2
#-------------------------------------------------------------
FGaxis<-rbind(faxis,gaxis)
FGcord <- rbind(frows, gcols)                                       # build a dataframe to be used as input for plotting via   
xmin <- min(FGcord[,firstaxis],FGcord[,lastaxis])
xmax <- max(FGcord[,firstaxis],FGcord[,lastaxis])
ymin <- min(FGcord[,lastaxis],FGcord[,firstaxis])
ymax <- max(FGcord[,lastaxis],FGcord[,firstaxis])
ellplot <- ggplot(FGcord, aes(x=FGcord[,firstaxis], y=FGcord[,lastaxis])) + 
  geom_point(aes(colour=categ, shape=categ), size=1.5) +
  geom_vline(xintercept = 0, linetype=2, color="gray") + 
  geom_hline(yintercept = 0, linetype=2, color="gray") + 
  labs(x=paste0("Dimension ", firstaxis,sep=" (", round(inertiapc[firstaxis],1),"%)"),  y=paste0("Dimension ", lastaxis,sep=" (", round(inertiapc[lastaxis],1),"%)" ))  +  
  scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin, ymax)) + 
  theme(panel.background = element_rect(fill="white", colour="black")) + 
  scale_color_manual(values=c("blue", "red")) + 
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE) + 
  geom_text_repel(data=FGcord, aes(colour=categ, label = labels), size = 3) +
  theme(legend.position="none")+
geom_ellipse(data=NULL,aes(x0=FGcord[,1],y0=FGcord[,2],a=FGaxis[,1],b=FGaxis[,2],angle=0,col=categ), stat = "ellip",
           position = "identity", n = 360, na.rm = FALSE, show.legend = NA,lty=2,
           inherit.aes = TRUE)+ 
ggtitle(" ") 
grid.arrange(ellplot, ncol=1)
#-------simple example  
#ggplot() +
#  geom_ellipse(aes(x0 = 0, y0 = 0, a = 1, b = 3, angle = 0, m1 = 3)) +
#  coord_fixed()
#-------------------------------------------------------
#xmin <- min(gcols[,firstaxis],gcols[,lastaxis])
#xmax <- max(gcols[,firstaxis],gcols[,lastaxis])
#ymin <- min(gcols[,lastaxis],gcols[,firstaxis])
#ymax <- max(gcols[,lastaxis],gcols[,firstaxis])
xmin <- min(gcols[,firstaxis],gcols[,lastaxis],hlax1.col,hlax2.col)
xmax <- max(gcols[,firstaxis],gcols[,lastaxis],hlax1.col,hlax2.col)
ymin <- min(gcols[,firstaxis],gcols[,lastaxis],hlax1.col,hlax2.col)
ymax <- max(gcols[,firstaxis],gcols[,lastaxis],hlax1.col,hlax2.col)
ellcol<-ggplot(gcols,aes(x=gcols[,firstaxis], y=gcols[,lastaxis]))+
  geom_point(aes(colour=categ, shape=categ), size=1.5) +
  geom_vline(xintercept = 0, linetype=2, color="gray") + 
  geom_hline(yintercept = 0, linetype=2, color="gray") + 
  labs(x=paste0("Dimension ", firstaxis,sep=" (", round(inertiapc[firstaxis],1),"%)"),  y=paste0("Dimension ", lastaxis,sep=" (", round(inertiapc[lastaxis],1),"%)" ))  +  
 # scale_x_continuous(limits = c(xmin, xmax)) +
#  scale_y_continuous(limits = c(ymin, ymax)) + 
   theme(panel.background = element_rect(fill="white", colour="black")) + 
  geom_text_repel(data=gcols, aes(colour=categ, label = labels), size = 3) +
  scale_color_manual(values=c("blue", "blue")) + 
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE) + 
        theme(legend.position="none")+     
    geom_ellipse(data=NULL,aes(x0=gcols[,1],y0=gcols[,2],a=hlax1.col,b=hlax2.col,angle=0,col="blue"), stat = "ellip",
             position = "identity", n = 360, na.rm = FALSE, show.legend = NA,
             inherit.aes = TRUE,lty=2)+
ggtitle(" ") 
grid.arrange(ellcol, ncol=1)

# ----------------------------------------------------ellipses on row categories
xmin <- min(frows[,firstaxis],frows[,lastaxis])
xmax <- max(frows[,firstaxis],frows[,lastaxis])
ymin <- min(frows[,lastaxis],frows[,firstaxis])
ymax <- max(frows[,lastaxis],frows[,firstaxis])
ellrow<-ggplot(frows,aes(x=frows[,firstaxis], y=frows[,lastaxis]))+
  geom_point(aes(colour=categ, shape=categ), size=1.5) +
  geom_vline(xintercept = 0, linetype=2, color="gray") + 
  geom_hline(yintercept = 0, linetype=2, color="gray") + 
  labs(x=paste0("Dimension ", firstaxis,sep=" (", round(inertiapc[firstaxis],1),"%)"),  y=paste0("Dimension ", lastaxis,sep=" (", round(inertiapc[lastaxis],1),"%)" ))  +  
  #scale_x_continuous(limits = c(xmin, xmax)) +
  #scale_y_continuous(limits = c(ymin, ymax)) + 
   theme(panel.background = element_rect(fill="white", colour="black")) + 
  geom_text_repel(data=frows, aes(colour=categ, label = labels), size = 3) +
  scale_color_manual(values=c("red", "red")) + 
    theme(legend.position="none")+  
  geom_ellipse(data=NULL,aes(x0=frows[,1],y0=frows[,2],a=hlax1.row,b=hlax2.row,angle=0,col="red"), stat = "ellip",
             position = "identity", n = 50, na.rm = FALSE, show.legend = FALSE,
             inherit.aes = TRUE,lty=2)+
ggtitle(" ") 
grid.arrange(ellrow, ncol=1)

#---------------------------------p-values
eccentricity <- abs(1 - (dmu[2]^2/dmu[1]^2))^(1/2)
area.row <- vector(mode = "numeric", length = I)
area.col <- vector(mode = "numeric", length = J)
for (i in 1:I) {
  area.row[i] <- 3.14159 * hlax1.row[i] * hlax2.row[i]
}
for (j in 1:J) {
  area.col[j] <- 3.14159 * hlax1.col[j] * hlax2.col[j]
}
pvalrow <- vector(mode = "numeric", length = I)
pvalcol <- vector(mode = "numeric", length = J)
for (i in 1:I) {
  if (M > 2) {
    pvalrow[i] <- 1- pchisq(t.inertia *n* ((1/Imass[i,i] - sum(a[i, 
                                                                 3:M]^2))^(-1)) * (cord1[i, 1]^2/dmu[1]^2 + cord1[i, 
                                                                                                               2]^2/dmu[2]^2), df = (I - 1) * (J - 1))
  }
  else {
    pvalrow[i] <- 1-pchisq(t.inertia*n * (Imass[i,i]) * 
                             (cord1[i, 1]^2/dmu[1]^2 + cord1[i, 2]^2/dmu[2]^2), 
                           df = (I - 1) * (J - 1))
  }
}
for (j in 1:J) {
  if (M > 2) {
    pvalcol[j] <-  1-pchisq(t.inertia *n* ((1/Jmass[j,j] - 
                                              sum(b[j, 3:M]^2))^(-1)) * (cord2[j, 1]^2/dmu[1]^2 + 
                                                                           cord2[j, 2]^2/dmu[2]^2), df = (I - 1) * (J - 1))
  }
  else {
    pvalcol[j] <- 1- pchisq(t.inertia *n* (Jmass[j,j]) * 
                              (cord2[j, 1]^2/dmu[1]^2 + cord2[j, 2]^2/dmu[2]^2), 
                            df = (I - 1) * (J - 1))
  }
}
summ.name <- c("HL Axis 1", "HL Axis 2", "Area", "P-value")
row.summ <- cbind(hlax1.row, hlax2.row, area.row, pvalrow)
dimnames(row.summ) <- list(paste(Inames), paste(summ.name))
col.summ <- cbind(hlax1.col, hlax2.col, area.col, pvalcol)
dimnames(col.summ) <- list(paste(Jnames), paste(summ.name))
(list(eccentricity = eccentricity, row.summ = row.summ, 
               col.summ = col.summ))
}