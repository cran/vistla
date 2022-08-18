#' Synthetic data representing a simple mediator chain
#'
#' Chain is generated from a simple Bayes network,
#' \deqn{X\rightarrow M_1 \rightarrow M_2  \rightarrow M_3  \rightarrow M_4  \rightarrow Y}
#' where every variable is binary.
#' The set consists of 11 observations, and is tuned to be easily deciphered.
#' @format A data set with six columns, each is a factor of two levels.
#' @usage data(chain)
"chain"

#' Synthetic data representing a junction
#'
#' Junction is a model of a multimodal agent, a variable that is an element of multiple separate paths.
#' Here, these paths are \eqn{A_1\rightarrow X \rightarrow A_2} and 
#' \eqn{B_1\rightarrow X \rightarrow B_2,}
#' while \eqn{X} is the junction.
#' The set consists of 12 observations, and is tuned to be easily deciphered.
#' @format A data set with five columns, each is a factor of two or four levels.
#' @usage data(junction)
"junction"

