#' @useDynLib vistla, .registration=TRUE
NULL

#' @rdname vistla
#' @export
vistla<-function(x,...)
 UseMethod("vistla")

#' @rdname vistla
#' @param formula alternatively, formula describing the task, in a form \code{root~predictors}, which adheres to standard R behaviours.
#'  Accepts \code{+} to add a predictor, \code{-} to omit one, and \code{.} to import whole \code{data}.
#'  Use \code{\link{I}} to calculate new predictors.
#'  When present in \code{data}, response is getting omitted from predictors.
#' @param data \code{data.frame} in context of which the formula will be executed; can be omitted when not using \code{.}.
#' @export
#' @method vistla formula
vistla.formula<-function(formula,data,...,yn){
 make_df(formula,data,parent.frame())->x
 if(missing(yn)) yn<-x$yn
 vistla.data.frame(x$X,x$Y,yn=yn,...)
}

#' Influence path identification with the Vistla algorithm
#'
#' Detects influence paths.
#' @rdname vistla
#' @param x data frame of predictors.
#' @param y vistla tree root, a feature from which influence paths will be traced.
#' @param flow algorithm mode, specifying the iota function which gives local score to an edge of an edge graph. If in doubt, use the default, \code{"fromdown"}.
#' @param iomin score threshold below which path is not considered further. The higher value the less paths are generated, which also lowers the time taken by the function.
#'  The default value of 0 turns of this filtering.
#'  The same effect can be later achieved with the \code{\link{prune}} function.
#' @param targets a vector of target feature names.
#'  If given, the algorithm will stop just after reaching the last of them, rather than after tracing all paths from the root.
#'  The same effect can be later achieved with the \code{\link{prune}} function.
#' @param verbose when set to \code{TRUE}, turns on reporting of the algorithm progress.
#' @param yn name of the root (\code{Y} value), used in result pretty-printing and plots. Must be a single-element character vector.
#' @param threads number of threads to use. 
#'  Value of 0 indicates all available for OpenMP.
#' @param ... pass-through arguments, ignored.
#' @return The tracing results represented as an object of a class \code{vistla}.
#'  Use \code{\link{paths}} and \code{\link{path_to}} functions to extract individual paths,
#'  \code{\link{branches}} to get the whole tree and \code{\link{mi_scores}} to get the basic score matrix.
#' @method vistla data.frame
#' @export
vistla.data.frame<-function(x,y,...,flow=c("fromdown","intoup","both","spread","from","into","up","down"),iomin=0,targets,verbose=FALSE,yn="Y",threads=0L){
 if(missing(targets)){
  targets<-integer(0)
 }else{
  targets<-unique(match(targets,names(x)))
  if(any(is.na(targets))) stop("Unknown variables specified as targets")
 }
 ans<-.Call(
  C_vistla,x,y,
  flow=switch(match.arg(flow),from=2L,into=1L,both=3L,spread=0L,fromdown=2L+8L,intoup=1L+4L,up=4L,down=8L),
  iomin,targets,verbose,
  as.integer(threads)
 )
 ans$yn<-yn
 ans$iomin<-iomin
 stats::setNames(
  data.frame(ans$tree),
  c("a","b","c","score","depth","leaf","used","prv")
 )->ans$tree

 if(length(targets)>0)
  ans<-prune_targets(ans,targets)

 class(ans)<-"vistla"
 ans
}

#' @rdname vistla
#' @method vistla default
#' @export
vistla.default<-function(x,...)
 stop("Expecting a formula or a data.frame as an input, got ",paste(class(x),collapse="/"))

#' Print vistla objects
#'
#' Utility functions to print vistla objects.
#' @method print vistla
#' @param x vistla object.
#' @param n maximal number of paths to preview.
#' @param ... ignored.
#' @return Invisible copy of \code{x}.
#' @export
print.vistla<-function(x,n=7L,...){
 stopifnot(inherits(x,"vistla"))
 stopifnot(n>0)
 pruned<-(x$iomin>0)||(!is.null(x$targets))
 cat(sprintf("\n\tVistla tree rooted in %s%s\n\n",x$yn,ifelse(pruned,", pruned","")))
 nl<-sum(x$tree$leaf)
 if(nl==0){
  cat("No paths.\n\n")
  return(invisible(x))
 }

 prev<-if(is.null(x$targets)){
  x$tree[x$tree$leaf & x$tree$used,c("c","score","depth")]
 }else{
  x$tree[(x$tree$c%in%x$targets) & x$tree$leaf & x$tree$used,c("c","score","depth")]
 }
 pn<-nrow(prev)
 if(pn>n) prev<-prev[1:n,]

 cat(
  if(is.null(x$targets)){
   if(pn==1) "Path:" else "Paths:"
  }else{
   if(pn==1) "Path to specified target:" else "Paths to specified targets:"
  }
 )

 prev$c<-colnames(x$mi)[prev$c]
 for(e in 1:nrow(prev)){
  pp<-path_to(x,prev$c[e])[-1]
  if(pn>1) pp[-length(pp)]->pp
  pp<-paste(pp,collapse=" ~ ")
  cat(paste(strwrap(sprintf(
   if(e==1) "%s (score %0.2g) ~ %s" else "%s (%0.2g) ~ %s",
   prev$c[e],prev$score[e],pp),initial='\n - ',prefix=' '),collapse="\n"))
 }
 if(pn>n) cat(sprintf("\n\t[%d more]",pn-nrow(prev)))
 if(pn>1) cat(sprintf("\n\t\t ... ~ %s\n\n",x$yn)) else cat("\n\n")
 invisible(x)
}

#' Coerce data as vistla would
#'
#' This function will coerce the input vector into factor as \code{\link{vistla}} function would.
#' Useful for testing or pre-computing quantisation.
#' @param x Input vector.
#' @return \code{x} coerced into factor.
#' @export
vistla_coerce<-function(x)
 .Call(C_convert,x)
