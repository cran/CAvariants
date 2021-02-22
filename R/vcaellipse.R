#library(ggplot2)
#library(ggforce)
#library(plotly)
#library(CAvariants)
vcaellipse<-function(row.summ,col.summ,inertiapc,cord1,cord2,a,b,firstaxis=1,lastaxis=2,eccentricity=eccentricity){
#globalVariables(categ)
thinglabels<-dimnames(cord1)[[1]]
varlabels<-dimnames(cord2)[[1]]
nthings<-dim(cord1)[1]
nvars<-dim(cord2)[1]
#---------------------------------------------------
eccentricity <- eccentricity
hlax1.row<-row.summ[,1]
hlax2.row<-row.summ[,2]
hlax1.col<-col.summ[,1]
hlax2.col<-col.summ[,2]
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
#-------------------------------------------------------
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

#---------------------------------results
(list(eccentricity = eccentricity, row.summ = row.summ, 
               col.summ = col.summ))
}