#' Construct the value for the flow argument
#'
#' Vistla builds the tree by optimising the influence score over path, which is given by the iota function.
#' The \code{flow} argument of the vistla function can be used to modify the default iota and some associated behaviours.
#' This function can be used to construct the proper value for this argument.
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
 codes<-c('from'=2L,'from!'=18L,'into'=1L,'into!'=17L,'spread'=0L,'spread!'=16L,
  'fromdown'=2L+8L,'both'=3L,'both!'=3L+16L,'intoup'=1L+4L,'down'=8L,'up'=4L)
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
   c(into,from,up,down,forcepath)*
   2^(0:4)
  ))
 }
 class(ans)<-"vistla_flow"
 ans
}

flow2char<-function(x)
 paste(
  c("spread","into","from","both")[bitwAnd(x,3L)+1],
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

