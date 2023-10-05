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
#' @param flow algorithm mode, specifying the iota function which gives local score to an edge of an edge graph. 
#'  If in doubt, use the default, \code{"fromdown"}.
#' @param iomin score threshold below which path is not considered further. 
#'  The higher value the less paths are generated, which also lowers the time taken by the function.
#'  The default value of 0 turns of this filtering.
#'  The same effect can be later achieved with the \code{\link{prune}} function.
#' @param targets a vector of target feature names.
#'  If given, the algorithm will stop just after reaching the last of them, rather than after tracing all paths from the root.
#'  The same effect can be later achieved with the \code{\link{prune}} function.
#'  This is a simple method to remove irrelevant paths, yet it comes with a substantial increase in computational burden.
#' @param verbose when set to \code{TRUE}, turns on reporting of the algorithm progress.
#' @param estimator mutual information estimator to use.
#'  \code{"mle"} --- maximal likelihood, requires all features to be discrete (factors or booleans).
#'  \code{"kt"} --- Kendall transformation, requires all features to be either ordinal (numeric, integer or ordered factor) or bi-valued (two-level factors or booleans).
#' @param yn name of the root (\code{Y} value), used in result pretty-printing and plots. 
#'  Must be a single-element character vector.
#' @param threads number of threads to use. 
#'  When missing or set to 0, vistla uses all available cores.
#' @param ... pass-through arguments, ignored.
#' @return The tracing results represented as an object of a class \code{vistla}.
#'  Use \code{\link{paths}} and \code{\link{path_to}} functions to extract individual paths,
#'  \code{\link{branches}} to get the whole tree and \code{\link{mi_scores}} to get the basic score matrix.
#' @references "Kendall transformation brings a robust categorical representation of ordinal data" M.B. Kursa. SciRep 12, 8341 (2022).
#' @method vistla data.frame
#' @export
vistla.data.frame<-function(x,y,...,flow,iomin,targets,estimator=c("mle","kt"),verbose=FALSE,yn="Y",threads){
 targets<-if(!missing(targets))
  unique(match(targets,names(x))) else integer(0)
 if(any(is.na(targets))) stop("Unknown variables specified as targets")
 if(missing(flow)) flow<-1L+8L
 if(is.character(flow)) flow<-vistla::flow(flow)
 if(missing(iomin)) iomin<-0
 if(missing(threads)) threads<-0L

 #Prepare input
 estimator<-match.arg(estimator)
 ec<-if(estimator=="mle") 1L else if(estimator=="kt") 2L else 17L;
 
 #Execute
 ans<-.Call(
  C_vistla,x,y,
  as.integer(flow)[1],ec,
  iomin,targets,verbose,
  as.integer(threads)
 )

 #Enhance the output
 ans$yn<-yn
 ans$iomin<-iomin
 stats::setNames(
  data.frame(ans$tree),
  c("a","b","c","score","depth","leaf","used","prv")
 )->ans$tree

 class(ans)<-"vistla"
 
 ans$flow<-flow
 class(ans$flow)<-"vistla_flow"
 
 if(length(targets)>0)
  ans<-prune_targets(ans,targets)
    
 ans
}

#' Construct the value for the flow
#'
#' Vistla builds the tree by optimising the influence score over path, which is given by the iota function.
#' The \code{flow} argument of the vistla function can be used to modify the default iota and some associated behaviours.
#' This function can be used to construct the proper value of this argument.
#' @param code Character code of the flow parameter, like \code{"fromdown"}. 
#'  If given, overrides other arguments.
#' @param from if \code{TRUE}, paths must satisfy data processing inequality as going from the root.
#' @param into if \code{TRUE}, paths must satisfy data processing inequality as going into the root.
#' @param down if \code{TRUE}, subsequent features on the path must have lower mutual information with the root; by default, true when \code{from} is true but if both \code{from} and \code{into} are true.
#' Can't be true together with \code{up}.
#' @param up if \code{TRUE}, subsequent features on the path must have higher mutual information with the root; by default, true when \code{into} is true but if both \code{from} and \code{into} are true.
#' Can't be true together with \code{down}.
#' @param forcepath when neither \code{up} or \code{down} is true, vistla may output walks rather than paths, i.e., sequences of features which are not unique.
#' Yet, when this argument is set to \code{TRUE}, additional condition is checked to forbid such self-intersections.
#' One should note that this check is computationally expensive, though.
#' By default true when both \code{up} and \code{down} are false.
#' @param ... ignored.
#' @return A \code{vistla_flow} object which can be passed to the \code{vistla} function; 
#'  in practice, a single integer value.
#' @export
flow<-function(code,...,from=TRUE,into=FALSE,down,up,forcepath){
 codes<-c('from'=1L,'from!'=17L,'into'=2L,'into!'=18L,'spread'=0L,'spread!'=16L,
  'fromdown'=1L+8L,'both'=3L,'both!'=3L+16L,'intoup'=2L+4L,'down'=8L,'up'=4L)
 if(!missing(code)){
  ans<-if(is.integer(code)) code else {
   if(is.character(code)){
    codes[code]
   }else stop("Wrong value of code ")
  }
  if(!is.integer(ans) || is.na(ans)) stop("Unknown code ",code)
 }else{
  if(missing(down)) down<-from & !into
  if(missing(up)) up<-into & !from
  if(missing(forcepath)) forcepath<-!up & !down
  ans<-as.integer(sum(
   c(from,into,up,down,forcepath)*
   2^(0:4)
  ))
 }
 class(ans)<-"vistla_flow"
 ans
}

flow2char<-function(x)
 paste(
  c("spread","from","into","both")[bitwAnd(x,3L)+1],
  ifelse(bitwAnd(x,4L)>0,"up",""),
  ifelse(bitwAnd(x,8L)>0,"down",""),
  ifelse(bitwAnd(x,16L)>0,"!",""),
  sep=""
 )

#' @rdname flow
#' @method print vistla_flow
#' @param x flow value to print.
#' @export
print.vistla_flow<-function(x,...){
 cat(paste(
  "Vistla flow: ",
  flow2char(x),
  "\n",
  sep=""))
 invisible(x)
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

