#' Measure agreement between vistla trees
#'
#' @param x a \code{vistla} object (first to compare) or a list of \code{vistla} objects (all to compare pairwise).
#' @param y a \code{vistla} object, second to compare, if \code{x} is a single object.
#' @param ... ignored.
#' @param method correlation method to use for quantification.
#'  See \code{\link{cor}} for possible values.
#' @param raw if \code{TRUE}, suppresses correlation calculation and output the raw aligned scores instead.
#' @return Correlation matrix with score correlations between each pair of given vistla trees.
#' @examples
#'  data(chain)
#'  agreement(
#'   vistla(Y~.,data=chain),
#'   vistla(Y~.,data=chain[,sample(6)])
#'  )
#' @export
agreement<-function(x,y=NULL,...,method="spearman",raw=FALSE){
 x<-if(is.null(y) && is.list(x)) x else list(x=x,y=y)
 stopifnot(all(sapply(x,inherits,"vistla")))
 if(length(x)==1) stop("Comparison requires at least two vistla objects")
 n<-names(x)
 KEY<-c('a','b','c')
 x<-lapply(lapply(x,data.frame),'[',,c(KEY,'score'))
 for(e in 1:length(x)) names(x[[e]])[4]<-sprintf("x%d",e)
 A<-NULL
 for(e in x)
  A<-if(is.null(A)) e else merge(A,e,by=KEY,all=TRUE,no.dups=FALSE)

 if(raw){
  if(length(n)==(length(names(A))-3))
   names(A)[-(1:3)]<-n
  return(A)
 }

 A<-as.matrix(A[,-(1:3)])
 A[is.na(A)]<-0 #Zero not-matched scores
 colnames(A)<-n
 if(ncol(A)==2) stats::cor(A[,1],A[,2],method=method) else stats::cor(A,method=method)
}
