#' Collapse the vistla tree into a pairwise graph
#'
#' @param x vistla object or a vistla_hierarchy object to collapse.
#' @param aggregate score aggregation mode.
#' "max" is the maximal score for this edge over all paths in the tree.
#' For raw vistla scores it means the score of the widest path this edge was a part of;
#' for ensemble scores, it corresponds to the count of the most often appearing path with this edge.
#' "sum" is the sum of scores. Makes little sense for raw vistla scores; for ensemble scores it corresponds
#' to the total count of this edge over all paths in the ensemble.
#' "none" returns a vector of scores over all paths, which can be processed anyhow the user desires.
#' @return A pairlist representation of the graph resulting from the tree collapse.
#'  The result is a data frame with the following columns.
#'  \code{A} & \code{B} are the ends of the edge, in order where A is closer to root than B
#'  (interpretation depends on the \code{flow} parameter used in \code{\link{vistla}} invocation);
#'  \code{score} is the score aggregated according to the \code{aggregate} argument;
#'  finally \code{paths} is the count of paths which included this edge.
#' @examples
#' \dontrun{
#'  data(junction)
#'  v<-vistla(Y~.,data=junction)
#'  collapse(v)
#' }
#' @export
collapse<-function(x,aggregate=c("max","sum","none")){
 if(inherits(x,"vistla")) x<-hierarchy(x)
 if(!inherits(x,"vistla_hierarchy")) stop("Collapse can only work on vistla or vistla_hierarchy objects")
 x$prv_name<-x$name[x$prv]
 x[x$depth>=0,]->x
 paste(
  as.numeric(factor(x$name)),
  as.numeric(factor(x$prv_name))
 )->pair_id
 aggregate<-match.arg(aggregate)
 
 AGG<-switch(aggregate,
  max=max,
  sum=sum,
  none=function(x) I(list(x))
 )
 do.call(rbind,lapply(split(x,pair_id),function(ch){
  data.frame(
   A=ch$prv_name[1],
   B=ch$name[1],
   score=AGG(ch$score),
   paths=nrow(ch)
  )
 }))->ans
 rownames(ans)<-NULL
 ans
}

