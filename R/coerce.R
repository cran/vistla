ew_cut<-function(bins)
 function(x) factor(paste('l',cut(x,breaks=bins,labels=FALSE),sep=""))

es_cut<-function(bins)
 function(x) ew_cut(bins)(rank(x,na.last=FALSE))


#' Basic discretisation of numerical features
#'
#' One can use this function for a quick, ad hoc discretisation of numerical features in a data frame, so that it could be passed to \code{\link{vistla}} using the maximal likelihood estimation (mle, the default).
#' This can be used to simulate legacy behaviour of vistla, which was to automatically perform such conversion with 10 equal-width bins.
#' The non-numeric columns are left as they were, hence this function is idempotent and does nothing when given fully discrete data.
#' @param x Data frame to be converted.
#' @param bins Number of bins to cut each numerical column into.
#' @param equal If given \code{"width"}, function performs cuts into bins of an equal width, which may thus contain substantially different number of objects.
#' One the other hand, when given \code{"size"} (default), cuts are done according to quantiles, hence provide bins with approximately the same number of objects, yet with different widths.
#' Both options are asymptotically equivalent when the distribution of a given column is uniform.
#' @return A copy of \code{x}, in which numerical columns have been discretised.
#' @note While convenient, this function does not necessary provide optimal quantisation of the data (in terms of future vistla performance); especially the bins parameter should be adjusted to the input data, either via optimisation or based on the known properties of the input or mechanisms behind it.
#' @examples
#' \dontrun{
#' data(cchain)
#' vistla(Y~.,data=mle_coerce(cchain,3,"size")) 
#' }
#' @export
mle_coerce<-function(x,bins=3,equal=c("size","width")){
 stopifnot(is.data.frame(x))
 equal<-match.arg(equal)
 qf<-if(equal=="width") ew_cut(bins) else 
      if(equal=="size") es_cut(bins) else stop("Unknown value of the `equal` argument")
 for(e in 1:ncol(x))
  if(is.numeric(x[,e])) x[,e]<-qf(x[,e])
 return(x)
}
