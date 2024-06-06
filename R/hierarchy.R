#' Extract the vertex hierarchy from the vistla tree
#'
#' Traverses the vistla tree in a depth-first order and 
#' lists the visited vertices as a data frame.
#' @param x vistla object.
#' @return A data frame of a class \code{vistla_hierarchy}.
#' @note This function effectively prunes the tree off suboptimal paths.
#' @export
hierarchy<-function(x){
 stopifnot(inherits(x,"vistla"))
 if(sum(x$tree$leaf)==0){
  Q<-data.frame(
   name=x$yn,
   score=NA_real_,
   depth=-1,
   leaf=FALSE,
   prv=NA_integer_
  )
  names(Q)<-c("name","score","depth","leaf","prv")
  class(Q)<-c("vistla_hierarchy","data.frame")
  return(Q)
 }
 Q<-subset_tree(x$tree,x$tree$used)
 which(is.na(Q$a))->ti
 Qt<-Q[ti,]
 Qt[Qt$used,]->Qt
 Qt[!duplicated(Qt$b),]->Qt
 nt<-nrow(Qt)
 Qy<-data.frame(a=NA,b=NA,c=NA,score=NA,depth=-1,leaf=FALSE,used=TRUE,prv=NA)
 Qh<-data.frame(
  a=rep(NA,nt),
  b=rep(NA,nt),
  c=Qt$b,
  score=Qt$score,
  depth=0,
  leaf=FALSE,
  used=TRUE,
  prv=1
 )
 Q$prv<-Q$prv+nt+1
 Q$prv[ti]<-match(Q$b[ti],Qh$c)+1
 rbind(Qy,Qh,Q)->Q

 Q$firstsub<-Q$lastsub<-Q$nxtsib<-rep(NA,nrow(Q))
 for(e in 1:nrow(Q)){
  prv<-Q$prv[e]
  if(!is.na(prv)){
   if(is.na(Q$firstsub[prv])) Q$firstsub[prv]<-e
   if(!is.na(Q$lastsub[prv])) Q$nxtsib[Q$lastsub[prv]]<-e
   Q$lastsub[prv]<-e
  }
 }

 Q$trord<-rep(NA,nrow(Q))
 Q$visited<-rep(FALSE,nrow(Q))
 cur<-1
 ord<-1
 while(!is.na(cur))
  cur<-if(Q$visited[cur]){
   if(!is.na(Q$nxtsib[cur])) Q$nxtsib[cur] else Q$prv[cur]
  }else{
   Q$trord[ord]<-cur
   Q$visited[cur]<-TRUE
   ord<-ord+1
   if(!is.na(Q$firstsub[cur])) 
    Q$firstsub[cur] else
     if(!is.na(Q$nxtsib[cur])) Q$nxtsib[cur] else Q$prv[cur]
  }
 
 Q<-subset_tree(Q,Q$trord)
 Q[,c("c","score","depth","leaf","prv")]->Q
 names(Q)[1]<-"name"
 Q$name<-colnames(x$mi)[Q$name]
 Q$name[is.na(Q$name)]<-x$yn
 rownames(Q)<-NULL
 class(Q)<-c("vistla_hierarchy","data.frame")

 return(Q)
}

#' @rdname print.vistla
#' @method print vistla_hierarchy
#' @export
print.vistla_hierarchy<-function(x,...){
 cat("\n\tVistla hierarchy\n\n")
 V<-'| '
 J<-'+'
 cat(sprintf(
   "%s%s%s%s",
   sapply(x$depth+1,function(dc) paste(rep(V,dc),collapse="")),
   J,x$name,
   ifelse(x$depth>0,sprintf(" (%0.4g)",x$score),"")
  ),sep="\n")
 cat("\n")
 invisible(x)
}

