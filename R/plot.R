#' Overview plot of the vistla tree
#'
#' @method plot vistla
#' @param x vistla object.
#' @param scale_width if \code{TRUE}, widths of links are scaled according to score.
#' @param ... additional graphical parameters, passed to plot.
#' @return \code{x}, invisibly.
#' @export
plot.vistla<-function(x,...,scale_width=TRUE){
 organise_plot(x)->trg
 yw<-graphics::strwidth(trg$name[trg$depth==-1],'inches')
 ew<-max(graphics::strwidth(trg$name[trg$depth==max(trg$depth)],'inches'))
 op<-graphics::par(omi=c(0,0,0,0),mai=c(0,yw/2,0,ew/2+0),xpd=NA)
 on.exit(graphics::par(op))
 plot(trg$x,trg$y,type="n",xlab="",yaxt='n',xaxt='n',ylab='',bty='n',...)
 trg$wd<-if(scale_width) 
  (trg$score-min(trg$score,na.rm=TRUE))/(max(trg$score,na.rm=TRUE)-min(trg$score,na.rm=TRUE))*4.5+0.5
 else 1
 trg$wd[!is.finite(trg$wd)]<-1
 graphics::segments(trg$x,trg$y,trg$xp,trg$yp,col="#cccccc",lty=1,lwd=trg$wd)
 graphics::text(trg$x,trg$y,ifelse(trg$leaf,trg$name,sprintf("(%s)",trg$name)))
 invisible(x)
}

organise_plot<-function(x){
 hierarchy(x)->tr

 #Traversing from leaves to Y to reserve size;
 # height(v)<=max(sum(height(subs(v))),1)
 tr$H<-rep(0,nrow(tr))
 for(e in nrow(tr):1){
  tr$H[e]<-max(1,tr$H[e])
  #Push sizes back to super
  prv<-tr$prv[e]
  tr$H[prv]<-tr$H[prv]+tr$H[e]
 }

 tr$x<-tr$depth
 tr$yc<-rep(0,nrow(tr))
 for(e in 1:nrow(tr)){
  tr$prv[e]->prv
  if(is.na(prv)){
   tr$yc[e]<-0
   tr$y[e]<-tr$H[e]/2
  }else{
   tr$yc[e]<-tr$yc[prv]
   tr$y[e]<-tr$yc[prv]+tr$H[e]/2
   tr$yc[prv]<-tr$yc[prv]+tr$H[e]
  }
 }
 tr$y<-max(tr$y)-tr$y

 tr$H<-NULL
 tr$yc<-NULL

 tr$yp<-tr$y[tr$prv]
 tr$xp<-tr$x[tr$prv]

 data.frame(tr)
}
