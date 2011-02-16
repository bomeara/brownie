\name{treeHeight-methods}
\docType{methods}
\alias{treeHeight}
\alias{treeHeight-methods}
\alias{treeHeight,phylo4-method}
\alias{treeHeight,phylo4d-method}
\alias{treeHeight,phylo4d_ext-method}
\title{ Get the height (depth) of a tree }
\description{
	Get the absolute height (or depth) of a given tree.  This function uses \code{tipdate.ci} to accomplish this by first converting the argument to class 'phylo' and then calling it with that converted value.
}
\section{Methods}{
\describe{

\item{\code{signature(x = "phylo4")}}{
Get the tree height
}

\item{\code{signature(x = "phylo4d")}}{
Get the tree height
}

\item{\code{signature(x = "phylo4d_ext")}}{
Get the tree height
}
}}
\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
