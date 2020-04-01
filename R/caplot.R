caplot<-function(frows,gcols,firstaxis,lastaxis,nseg,inertiapc,thingseg,col1,col2,col3,size1,size2){
#biplots for nominal variables
categ<-NULL
#frows <- data.frame(coord=cord1, labels=thinglabels, categ=rep("rows", nthings)) # build a dataframe to be used as input for plotting via ggplot2
#gcols <- data.frame(coord=cord2, labels=varlabels, categ=rep("cols", nvars)) # build a dataframe to be used as input for plotting via ggplot2
FGcord <- rbind(frows, gcols)                                       # build a dataframe to be used as input for plotting via   
xmin <- min(FGcord[,firstaxis],FGcord[,lastaxis])
xmax <- max(FGcord[,firstaxis],FGcord[,lastaxis])
ymin <- min(FGcord[,lastaxis],FGcord[,firstaxis])
ymax <- max(FGcord[,lastaxis],FGcord[,firstaxis])
 CAplot <- ggplot(FGcord, aes(x=FGcord[,firstaxis], y=FGcord[,lastaxis]),type="b") + 
    geom_point(aes(color=categ, shape=categ), size=size1) +
    geom_vline(xintercept = 0, linetype=2, color="gray") + 
    geom_hline(yintercept = 0, linetype=2, color="gray") + 
    labs(x=paste0("Dimension ",firstaxis,sep=" (", round(inertiapc[1],1), "%)"),y=paste0("Dimension ",lastaxis,sep= " (", round(inertiapc[2],1),"%)"))  +  
    scale_x_continuous(limits = c(xmin, xmax)) +
    scale_y_continuous(limits = c(ymin, ymax)) + 
theme(panel.background = element_rect(fill="white", colour="black")) + 
  #  scale_colour_manual(values=c("red", "blue")) + 
 scale_colour_manual(values=c(col1, col2)) +   
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE) + 
    geom_text_repel(data=FGcord, aes(colour=categ, label = labels), size = size2) +
  geom_segment(data=thingseg,aes(x=rep(0,c(nseg)),y=rep(0,c(nseg)),xend=thingseg[, firstaxis],yend=thingseg[, 
lastaxis]),colour=rep(col3,nseg))+
theme(legend.position="none")+
   ggtitle(" ") 
  grid.arrange(CAplot, ncol=1)
 }
