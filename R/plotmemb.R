plotmemb <- function(X,cl,inorder=FALSE,main="",drawborder=TRUE,lwd=0.5,lwdb=1,lcol=2)
{
 require(lattice)
 require(Matrix)
 if (inorder) {
     tmp <- sort(table(cl),dec=TRUE)
 } else {
     tmp <- table(cl)
 }
 tmp <- as.numeric(names(tmp))
 ind2 <- unlist(sapply(1:max(tmp), function(x) which(cl==tmp[x])))
 X2 <- as.matrix(X)
 AS <- Matrix(X2[ind2,ind2])
 p <- image(AS,lwd=lwd,main=main)
 print(p)

 if (drawborder) {
     trellis.focus("panel", 1, 1, highlight=FALSE)
     if (length(lcol)==1) lcol <- rep(lcol,max(cl))
     nc <- table(cl)
     nc.x <- nc.y <- c(1,cumsum(nc)+1)
     for (i in 1:max(cl)) {
         x1 <- nc.x[i]; x2 <- nc.x[i+1]
         y1 <- nc.y[i]; y2 <- nc.y[i+1]
         lpolygon(c(x1,x1,x2,x2), c(y1,y2,y2,y1),border=lcol[i],lwd=lwdb)
     }
     trellis.unfocus()
 }
}
