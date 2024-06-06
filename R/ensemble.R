#Ensemble-mode vistal, delegated from vistla.data.frame
vistla_ensemble<-function(x,y,flow,iomin,targets,ec,yn="Y",ensemble,threads){
 #Input was coerced upstream, except for ensemble

 #Short-cut, allow just a number to be passed as n
 if(is.numeric(ensemble) && length(ensemble)==1)
  ensemble<-ensemble(n=ensemble)
 
 ens<-with(ensemble,
  c(
   as.integer(ensemble$n),
   if(is.numeric(resample)) 
    pmax(3,pmin(round(resample*nrow(x)),n-1)) else
    ifelse(resample,0,nrow(x)),
   as.integer(ensemble$prune)
  )
 )
  
 #Execute
 ans<-.Call(
  C_vistlaEnsemble,x,y,
  as.integer(flow)[1],ec,
  iomin,targets,
  as.integer(ens),
  as.integer(threads)
 )
 ans<-data.frame(stats::setNames(ans,c("name","score","depth","leaf","prv")))
 ans$name<-names(x)[ans$name]
 attr(ans,"ensemble")<-ensemble
 ans$name[1]<-yn
 class(ans)<-c("vistla_hierarchy","data.frame")

 if(length(targets)>0) ans<-prune(ans,targets=names(x)[targets])

 ans
}

#' Construct the value for the ensemble argument
#'
#' Vistla can be run in the ensemble mode, in which tree is built multiple times, usually on a slightly modified input data.
#' This mode can be triggered by passing a value to the ensemble argument of the vistla method.
#' This function can be used to construct the proper value for this argument.
#' @param n number of replicatons.
#' @param resample if \code{TRUE}, a modified bootstrap is used; that is, algorithm draws as many objects as are in the original data, but with replacement, hence only about 63.2% of objects remain. Though, in vistla objects are deduplicated to avoid estimation quirks.
#' If this argument is given a number, it is interpreted to randomly sample exacly this fraction of objects, without replacement.
#' Fraction \code{f} of \code{n} objects is interpreted as \code{round(n*f)}, but not less than 3 and no more than \code{n-1}.
#' If \code{FALSE}, no resampling is done (vistla trees are just built using different random seeds.
#' @param prune Minimal number of iterations in which certain branch must appear not be prunned during ensemble consolidation.
#' Zero (default) means no prunning.
#' Note that \code{iomin} and \code{targets} arguments of the base algorithm can also be used to control the size of the resulting consensus tree.
#' @param ... ignored.
#' @return A \code{vistla_ensemble_control} object which can be passed to the \code{vistla} function.
#' @export
ensemble<-function(n=30,resample=TRUE,prune=0){
 bad_resample<-c(  
 length(resample)!=1,
  (is.numeric(resample) && (!is.finite(resample) || resample<=0 || resample>=1)),
  (is.logical(resample) && is.na(resample)),
  (!is.numeric(resample) && !is.logical(resample)))
 if(sum(bad_resample)>0) stop("Value of resample should be a single value, logical or numerical in (0;1) range")

 n<-as.integer(n)
 if((length(n)!=1) || n<1 || is.na(n)) stop("Value of n should be a single, positive integer")

 prune<-as.integer(prune)
 if((length(prune)!=1) || prune<0 || prune>n) stop("Value of prune should be a single number no larger than n")
 ans<-list(n=n,resample=resample,prune=prune)
 class(ans)<-"vistla_ensemble_control"
 ans
}

#' @rdname ensemble
#' @method print vistla_ensemble_control
#' @param x ensemble control value to print.
#' @export
print.vistla_ensemble_control<-function(x,...){
 cat(
  sprintf(
   "Vistla ensemble control, %d replications%s.\n%s\n",
   x$n,
   ifelse(is.numeric(x$resample),
    sprintf(" with %0.6g%% subsampling",x$resample*100),
    ifelse(x$resample," with bootstrap","")
   ),
   if(!is.null(x$prune) && x$prune>0) sprintf("Consensus prune threshold: %d.\n",x$prune) else ""
  )
 )
 invisible(x)
}
