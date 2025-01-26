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
#'  Consult the documentation of the \code{\link{flow}} function for more information.
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
#' @param ensemble used to switch vistla to the ensemble mode, in which a number of vistla models are built over permuted realisations of the input, and merged into a single consensus tree.
#' Should be given an output of the \code{\link{ensemble}} function; as a short-cut, one can pass a single number, which will be interpreted as the number of replications with other ensemble parameter default.
#' That is, \code{ensemble=30} is equivalent to \code{ensemble=ensemble(n=30)}.
#' Permutations are applied before estimators.
#' @param ... pass-through arguments, ignored.
#' @return Normally, the tracing results represented as an object of a class \code{vistla}.
#'  Use \code{\link{paths}} and \code{\link{path_to}} functions to extract individual paths,
#'  \code{\link{branches}} to get the whole tree and \code{\link{mi_scores}} to get the basic score matrix.
#' 
#'  When \code{ensemble} argument is given, a hierarchy object with the scored being counts of times certain path was present among the replicated ensemble, possibly pruned. 
#' @note The ensemble mode is both faster and makes better use of multithreading than replicating vistla manually.
#' @references "Vistla: identifying influence paths with information theory" M.B. Kursa. Bioinformatics btaf036 (2025).
#' @references "Kendall transformation brings a robust categorical representation of ordinal data" M.B. Kursa. SciRep 12, 8341 (2022).
#' @method vistla data.frame
#' @export
vistla.data.frame<-function(x,y,...,flow,iomin,targets,estimator=c("mle","kt"),verbose=FALSE,yn="Y",ensemble,threads){
 targets<-if(!missing(targets))
  unique(match(targets,names(x))) else integer(0)
 if(any(is.na(targets))) stop("Unknown variables specified as targets")
 if(missing(flow)) flow<-10L
 if(is.character(flow)) flow<-vistla::flow(flow)
 if(missing(iomin)) iomin<-0
 if(missing(threads)) threads<-0L

 #Prepare input
 estimator<-match.arg(estimator)
 ec<-if(estimator=="mle") 1L else if(estimator=="kt") 2L else 17L;
 
 #Redirect flow in case ensembling is requested
 if(!missing(ensemble)) return(vistla_ensemble(x,y,flow,iomin,targets,ec,yn,ensemble,threads))
 
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

#' Extract mutual information score matrix
#'
#' Produces a matrix \eqn{S} where \eqn{S_{ij}} is a 
#'  value of \eqn{I(X_i;X_j)}.
#' This matrix is always calculated as an initial step of the
#'  vistla algorithm and stored in the vistla object.
#' @param x vistla object.
#' @return A symmetric square matrix with mutual information scores between features and root.
#' @export
mi_scores<-function(x){
 stopifnot(inherits(x,"vistla"))
 rbind(cbind(x$mi,x$miY),c(x$miY,NA))->ans
 colnames(ans)[length(x$miY)+1]<-x$yn
 rownames(ans)[length(x$miY)+1]<-x$yn
 ans
}

#' Extract leaf scores of vertex pairs
#'
#' Produces a matrix \eqn{S} where \eqn{S_{ij}} is a score
#'  of the path ending in vertices \eqn{i} and \eqn{j}.
#' Since vistla works on vertex pairs, this value is unique.
#' This can be interpreted as a feature similarity matrix
#'  in context of the current vistla root.
#' @note This function should be called on an unpruned vistla tree,
#'  otherwise the result will be mostly composed of zeroes.
#' @param x vistla object.
#' @return A square matrix with leaf scores of all feature pairs.
#' @export
leaf_scores<-function(x){
 stopifnot(inherits(x,"vistla"))
 x$mi*0->ans
 ans[as.matrix(x$tree[,c("b","c")])]<-x$tree$score
 ans
}
