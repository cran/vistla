% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/export.R
\name{write.dot}
\alias{write.dot}
\title{Export tree to a Graphviz DOT format}
\usage{
write.dot(
  x,
  con,
  vstyle = list(shape = function(x) ifelse(x$depth < 0, "egg", ifelse(x$leaf, "box",
    "ellipse")), label = function(x) sprintf("\\"\%s\\"", x$name)),
  estyle = list(penwidth = function(x) sprintf("\%0.3f", 0.5 + x$score/max(x$score) *
    2.5)),
  gstyle = list(overlap = "\\"prism\\"", splines = "true"),
  direction = c("none", "fromY", "intoY")
)
}
\arguments{
\item{x}{vistla object.}

\item{con}{connection; passed to \code{writeLines}.
If missing, the DOT code is returned as a character vector.}

\item{vstyle}{vertex attribute list --- should be a named list of Graphviz attributes like \code{shape} or \code{penwidth}. 
For elements which are strings or numbers, the value is copied as is as an attribute value.
For elements which functions, though, the function is called on a \code{vistla_tree} object and should return a vector of values.}

\item{estyle}{edge attribute list, behaves exactly like \code{vstyle}. 
When functions are called, the Y-vertex is not present.}

\item{gstyle}{graph attribute list. Functions are not supported here.}

\item{direction}{when set to \code{"none"}, graph is undirected, otherwise directed, for \code{"fromY"}, root is a source, while for \code{"intoY"}, a sink.}
}
\value{
For a missing \code{con} argument, a character vector with the graph in the DOT format, invisible \code{NULL} otherwise.
}
\description{
Exports the vistla tree in a DOT format, which can be later layouted and rendered by Graphviz programs like dot or neato.
}
\note{
Graphviz attribute values can be either strings, like \code{"some vertex"} in \code{label}, or atoms, like \code{box} for \code{shape}.
 When returning a string value, you must supply quotes, otherwise it will be included as an atom.

 The default value of \code{gstyle} may invoke long layout calculations in Graphviz.
 Change to \code{list()} for a fast but less aesthetic layout.

 The function does no validation whether provided attributes or values are correct.
}
\references{
"An open graph visualization system and its applications to software engineering" E.R. Gansner, S.C. North. Software: Practice and Experience 30:1203-1233 (2000).
}
