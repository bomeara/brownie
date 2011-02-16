\name{hasSubNodes-methods}
\docType{methods}
\alias{hasSubNodes}
\alias{hasSubNodes-methods}
\alias{hasSubNodes,list-method}
\alias{hasSubNodes,phylo-method}
\alias{hasSubNodes,phylo4-method}
\alias{hasSubNodes,phylo4d_ext-method}
\title{ Does an object contain subnodes? }
\description{
Check to see if a 'phylo4d_ext' objects actually contains subnodes
}
\section{Methods}{
\describe{

\item{\code{signature(x = "list")}}{
Apply this function to each object in the list
}

\item{\code{signature(x = "phylo")}}{
Returns FALSE
}

\item{\code{signature(x = "phylo4")}}{
Returns FALSE
}

\item{\code{signature(x = "phylo4d_ext")}}{
Checks if \code{subnode.data} contains any columns
}
}}
\keyword{methods}

