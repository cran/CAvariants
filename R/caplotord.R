caplotord<-function(frows,gcols,firstaxis,lastaxis,nseg,inertiapc,thingseg,col1,col2,col3,size1,size2){
#plotting biplot for ordered variables-depicting polynomials
#  FGcord <- rbind(frows, gcols)
categ<-NULL
  FGcord <- rbind(gcols,frows)
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
 scale_colour_manual(values=c(col2, col1)) +   
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE) + 
    geom_text_repel(data=FGcord, aes(colour=categ, label = labels), size = size2) +
geom_line(aes(group=categ,linetype=categ),lwd=.2,colour="gray")+
scale_linetype_manual(values=c("rows"="dashed","cols"="blank"))+
  geom_segment(data=thingseg,aes(x=rep(0,c(nseg)),y=rep(0,c(nseg)),xend=thingseg[, firstaxis],yend=thingseg[, 
lastaxis]),colour=rep(col3,nseg))+
theme(legend.position="none")+
   ggtitle(" ") 
  grid.arrange(CAplot, ncol=1)
 }
