#' Prune the vistla tree
#'
#' This function allows to filter out suboptimal branches, as well as weak ones or these not in particular paths of interest.
#' @param x vistla object or a vistla_hierarchy object.
#' @param targets a character vector of features.
#'  When not missing, all branches not on lying paths to these targets are pruned.
#'  Unreachable targets are ignored, while names not present in the analysed set cause an error.
#' @param score a score threshold below which branches should be removed.
#'  When given, it effectively overrides the value of \code{iomin} or \code{ensemble(prune,...)} given to the \code{vistla} invocation; to this end, it can only be higher then the original value, since prune only modifies the output and cannot re-run the pathfinding.
#' @param iomin a legacy name for score, valid only for vistla objects; passing a value to either of them works the same, but giving some values for both is an error.
#' @return Pruned \code{x}; if both arguments are missing, this function still removes suboptimal branches.
#' @examples
#' \dontrun{
#'  data(chain)
#'  v<-vistla(Y~.,data=chain)
#'  print(v)
#'  print(prune(v,targets="M3"))
#'  print(prune(v,score=0.3))
#' }
#' @export
prune<-function(x,targets,iomin,score){
 if(inherits(x,"vistla_hierarchy")){
  if(missing(targets)) targets<-unique(x$name[x$leaf])
  if(missing(score)) score<--Inf
  return(prune_hierarchy(x,targets,score))
 }
 if(!inherits(x,"vistla")) stop("Prune can only work on vistla or vistla_hierarchy objects")
 if(!missing(iomin)||!missing(score)){
  if(!missing(score))
   if(missing(iomin)) iomin<-score else stop("Score and iomin arguments cannot be given at the same time")
  if(iomin<x$iomin) stop("Prune can only increase iomin")
  subset_tree(x$tree,x$tree$score>iomin)->x$tree
  x$iomin<-iomin
 }
 if(!missing(targets)){
  match(targets,colnames(x$mi))->ti
  if(any(is.na(ti))) stop("Unknown names in targets")
  if(!is.null(x$targets)){
   #Tree was already pruned; are we trying to extend?
   if(!all((ti%in%x$targets)|(ti%in%(x$tree$c[x$tree$leaf&x$tree$used])))) stop("Prune can only remove targets")
  }
  return(prune_targets(x,ti))
 }
 prune_targets(x)
}

prune_hierarchy<-function(x,targets,score){
 ss<-x$score
 ss[is.na(ss)]<-Inf
 to_keep<-which(
  (x$leaf&(x$name%in%targets)&(ss>score))
  |(x$depth<0))
 keep<-rep(FALSE,nrow(x))
 while(length(to_keep)>0){
  keep[to_keep]<-TRUE
  x$prv[to_keep]->to_keep
  unique(to_keep[!is.na(to_keep)])->to_keep
 }
 subset_tree(x,keep)->ans
 rownames(ans)<-NULL
 ans
}

#Internal engine of prune, used by vistla.data.frame as well
prune_targets<-function(x,ti){
 tree<-x$tree
 to_keep<-if(missing(ti)) 
  which(tree$leaf) else which((tree$c%in%ti)&tree$leaf)
 keep<-rep(FALSE,nrow(tree))
 while(length(to_keep)>0){
  keep[to_keep]<-TRUE
  tree$prv[to_keep]->to_keep
  unique(to_keep[!is.na(to_keep)])->to_keep
 }
 x$tree<-subset_tree(tree,keep)
 if(nrow(x$tree)>0)
  x$tree$used<-TRUE
 if(!missing(ti)) x$targets<-ti
 x
}
