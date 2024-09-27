#' Synthetic data representing a simple mediator chain
#'
#' Chain is generated from a simple Bayes network,
#' \deqn{Y\rightarrow M_1 \rightarrow M_2  \rightarrow M_3  \rightarrow M_4  \rightarrow T}
#' where every variable is binary.
#' The set consists of 11 observations, and is tuned to be easily deciphered.
#' @format A data set with six binary factor columns.
#' @usage data(chain)
"chain"

#' Synthetic continuous data representing a simple mediator chain
#'
#' Chain is generated from an uniform variable X by progressively adding gaussian noise, producing a mediator chain identical to this of the \code{\link{chain}} data, i.e.,
#' \deqn{Y\rightarrow M_1 \rightarrow M_2  \rightarrow M_3  \rightarrow M_4  \rightarrow T}
#' The set consists of 20 observations, and is tuned to be easily deciphered.
#' @format A data set with six numerical columns.
#' @usage data(cchain)
"cchain"

#' Synthetic data representing a junction
#'
#' Junction is a model of a multimodal agent, a variable that is an element of multiple separate paths.
#' Here, these paths are 
#' \eqn{Y\rightarrow A_1\rightarrow A_2\rightarrow J \rightarrow A_3}
#' and 
#' \eqn{Y\rightarrow B_1\rightarrow B_2\rightarrow J \rightarrow B_3,}
#' while \eqn{J} is the junction.
#' The set consists of 50 observations.
#' @format A data set with eight factor columns.
#' @usage data(junction)
"junction"

