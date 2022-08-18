#' Extract all branches of the Vistla tree
#'
#' Gives access to a list of all branches in the tree.
#' @param x vistla object.
#' @param suboptimal if TRUE, sub-optimal branches are included.
#' @note Pruned trees (obatined with \code{\link{prune}} and using \code{targets} argument
#'  in the \code{\link{vistla}} call) have no suboptimal branches.
#' @return A data frame collecting all branches traced by vistla.
#'  Each row corresponds to a single branch, i.e., edge between feature pairs.
#'  This way it is a triplet of original features, names of which are stored in \code{a},
#'  \code{b} and \code{c} columns.
#'  For instance, path \eqn{I \rightarrow J \rightarrow K \rightarrow L \rightarrow M} 
#'  would be stored in three rows, for \eqn{(a,b,c)}=\eqn{(I,J,K)}, \eqn{(J,K,L)}
#'  and \eqn{(K,L,M)}.
#'  The width of a path (minimal \eqn{\iota} value) between root and feature pair \eqn{(b,c)} is
#'  stored in the \code{score} column.
#'  \code{depth} stores the path depth, starting from 1 for pairs directly connected to the root,
#'  and increasing by one for each additional feature.
#'  Final column, \code{leaf}, is a logical path indicating whether the edge is a final segment
#'  of the widest path between root and \eqn{c}.
#' @export
branches<-function(x,suboptimal=FALSE){
 stopifnot(inherits(x,"vistla"))
 tree<-x$tree[x$tree$used|suboptimal,c('a','b','c','score','depth','leaf')]
 tree$a<-colnames(x$mi)[tree$a]
 tree$a[is.na(tree$a)]<-x$yn
 tree$b<-colnames(x$mi)[tree$b]
 tree$c<-colnames(x$mi)[tree$c]

 rownames(tree)<-NULL
 tree
}

#' @rdname branches
#' @method as.data.frame vistla
#' @param row.names passed to \code{as.data.frame}.
#' @param optional passed to \code{as.data.frame}.
#' @param ... ignored.
#' @export
as.data.frame.vistla<-function(x,row.names=NULL,optional=FALSE,suboptimal=FALSE,...)
 as.data.frame(branches(x,suboptimal),row.names=row.names,optional=optional)

#' Extract a single path
#'
#' Gives access to a vector of feature names over a path to a certain target feature.
#' @param x vistla object.
#' @param detailed if \code{TRUE}, suppresses default output and presents the same paths as a data frame featuring score.
#' @param target target feature name.
#' @return By default, a character vector with names of features along the path from \code{target} into root.
#'  When \code{detailed} is set to \code{TRUE}, a \code{data.frame} in a format identical to this produced by
#'  \code{\link{branches}}, yet without the \code{leaf} column.
#' @export
path_to<-function(x,target,detailed=FALSE){
 stopifnot(inherits(x,"vistla"))
 colnames(x$mi)->n
 which(n==target)->idx
 if(length(idx)!=1) stop(sprintf("Target feature %s was not in the training set",target))
 which((x$tree$c==idx) & (x$tree$leaf))->br
 if(length(br)>0){
  ans<-integer(x$tree$depth[br])
  for(e in 1:length(ans)){
   ans[e]<-br
   br<-x$tree$prv[br]
  }
 }else{
  ans<-integer()
 }
 if(detailed){
  x$tree[ans,c('a','b','c','score')]->ans
  ans$a<-n[ans$a]
  ans$a[is.na(ans$a)]<-x$yn
  ans$b<-n[ans$b]
  ans$c<-n[ans$c]
  rownames(ans)<-NULL
 }else{
  if(length(ans)==0) return(character())
  n[c(x$tree$c[ans[1]],x$tree$b[ans],x$tree$a[ans[length(ans)]])]->ans
  ans[is.na(ans)]<-x$yn
 }
 ans
}

#' List all paths
#'
#' Executes \code{\link{path_to}} for all path possible targets and returns
#' a list with the results.
#' @param x vistla object.
#' @param targets_only if \code{TRUE}, only paths to targets are extracted. 
#'  By default, turned on when \code{x} has targets, and off otherwise.
#' @param detailed passed to \code{\link{path_to}}. If \code{TRUE},
#'  suppresses default output and presents the same paths in a form of
#'  data frames featuring score.
#' @return A named list with one element per leaf or target, containing
#'  the path between this feature and root, in a format identical
#'  to this used by the \code{\link{path_to}} function.
#' @export
paths<-function(x,targets_only=!is.null(x$targets),detailed=FALSE){
 targets<-if(targets_only) c(x$targets) else x$tree$c[x$tree$leaf&x$tree$used]
 colnames(x$mi)->n
 targets<-n[targets]
 stats::setNames(lapply(targets,function(t) path_to(x,t,detailed=detailed)),targets)
}

