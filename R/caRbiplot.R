caRbiplot <- function(frows, gcols, firstaxis, lastaxis, inertiapc,  bip="row", size1,size2){
  ##########################################################################
  # #
  # Principal and standard Coordinates # frows and gcols
  # #
  # #
  ##########################################################################
I<-nrow(frows)
J<-nrow(gcols) 
categ<-NULL
attnam<-NULL
slp<-NULL
#frows <- data.frame(coord=f, labels=dimnames(f)[[1]], categ=rep("rows", I), linet=rep("solid",I)) # build a dataframe to be used as input for plotting via ggplot2
#gcols <- data.frame(coord=g, labels=dimnames(g)[[1]], categ=rep("cols", J), linet=rep("blank",J)) # build a dataframe to be used as input for plotting via ggplot2
if (bip == "row"){
rglines=data.frame(d1=gcols[,firstaxis],d2=gcols[,lastaxis],attnam=gcols$categ)
#rglines=data.frame(d1=gcols[,firstaxis],d2=gcols[,lastaxis],attnam="blue")
rglines$slp=rglines$d2/rglines$d1
}
if ((bip == "column")||(bip == "col")) {
rglines=data.frame(d1=frows[,firstaxis],d2=frows[,lastaxis],attnam=frows$categ)
#rglines=data.frame(d1=frows[,firstaxis],d2=frows[,lastaxis],attnam="blue")
rglines$slp=rglines$d2/rglines$d1
}
#-------------------------------------------------------------
ndim1<-I
catall <- rep("solid", ndim1)
FGcord <- rbind(frows, gcols)                                       # build a dataframe to be used as input for plotting via   
xmin <- min(FGcord[,firstaxis],FGcord[,lastaxis])
 xmax <- max(FGcord[,firstaxis],FGcord[,lastaxis])
ymin <- min(FGcord[,lastaxis],FGcord[,firstaxis])
ymax <- max(FGcord[,lastaxis],FGcord[,firstaxis])
 CAplot <- ggplot(FGcord, aes(x=FGcord[,firstaxis], y=FGcord[,lastaxis]), type="b") + 
   geom_point(aes(color=categ), size=size1) +
scale_shape_manual(values=categ) +
    geom_vline(xintercept = 0, linetype=2, color="gray") + 
    geom_hline(yintercept = 0, linetype=2, color="gray") + 
    labs(x=paste0("Dimension", firstaxis,sep=" ", round(inertiapc[firstaxis],1), "%"),y=paste0("Dimension",lastaxis,sep= " ", round(inertiapc[lastaxis],1),"%"))  +  
    scale_x_continuous(limits = c(xmin, xmax)) +
    scale_y_continuous(limits = c(ymin, ymax)) + 
theme(panel.background = element_rect(fill="white", colour="black")) + 
   coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE) + 
    geom_text_repel(data=FGcord, aes(colour=categ, label = labels), size = size2, max.overlaps =Inf) +
theme(legend.position="none")+
#geom_abline(data=rglines,aes(intercept=0,slope=slp,colour="blue"),alpha=.5)
geom_abline(data=rglines,aes(intercept=0,slope=slp,colour=attnam),alpha=.5)
grid.arrange(CAplot, ncol=1)
#list( frows=frows, gcols = gcols)
}



