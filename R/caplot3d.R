#open3d()
caplot3d <- function(f, g, percIn, size=10, adj=c(0,0.5) ,prop=1){
picsize<-c(range(f[,c(1,2,3)], g[,c(1,2,3)])/prop)
#browser()
# plot3d(0, 0, 0, type = "n", box = FALSE, xlim = range(g[,1], g[,2], g[,3]),
#       ylim = range(g[, 1], g[, 2], g[, 3]),
#       zlim = range(g[, 1], g[, 2], g[, 3]),
plot3d(0, 0, 0, type = "n", box = FALSE, xlim =picsize, ylim=picsize,zlim=picsize, 
      xlab = paste("Principal Axis 1", "(", round(percIn[1], digits = 2), "%", ")"),
       ylab = paste("Principal Axis 2", "(", round(percIn[2], digits = 2), "%", ")"),
       zlab = paste("Principal Axis 3", "(", round(percIn[3], digits = 2), "%", ")"))
points3d(f[, 1], f[, 2], f[, 3], pch="*", size = size, col = "blue")
text3d(f[, 1], f[, 2], f[, 3], texts = dimnames(f)[[1]], adj = adj, col = "blue")
points3d(g[, 1], g[, 2], g[, 3], size = size, col = "red")
text3d(g[, 1], g[, 2], g[, 3], texts = dimnames(g)[[1]], adj = adj, col = "red")
i <- c(1, 2, 1, 3, 1, 4)
x1 <- c(0, max(abs(f), abs(g)), 0, 0)
y1 <- c(0, 0, max(abs(f),abs(g)), 0)
z1 <- c(0, 0, 0, max(abs(f), abs(g)))
x2 <- c(0, - max(abs(f), abs(g)), 0, 0)
y2 <- c(0, 0, - max(abs(f), abs(g)), 0)
z2 <- c(0, 0, 0, - max(abs(f), abs(g)))
segments3d(x1[i], y1[i], z1[i], alpha=0.3)
segments3d(x2[i], y2[i], z2[i], alpha=0.3)
#for (j in 1:6){
#  lines3d(c(0, g[j,1]), c(0,g[j,2]), c(0,g[j,3]), col="red", lty = 3)
#}
#rgl.postscript("foodcaplot3d.eps", fmt = "eps")
# Instead of using this command to make an eps file of the plot, save it into word,
#then convert to pdf. From pdf convert to an eps file
 }