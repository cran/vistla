% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{cchain}
\alias{cchain}
\title{Synthetic continuous data representing a simple mediator chain}
\format{
A data set with six numerical columns.
}
\usage{
data(cchain)
}
\description{
Chain is generated from an uniform variable X by progressively adding gaussian noise, producing a mediator chain identical to this of the \code{\link{chain}} data, i.e.,
\deqn{Y\rightarrow M_1 \rightarrow M_2  \rightarrow M_3  \rightarrow M_4  \rightarrow T}
The set consists of 20 observations, and is tuned to be easily deciphered.
}
\keyword{datasets}
