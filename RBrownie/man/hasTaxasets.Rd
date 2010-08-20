\name{hasTaxasets-methods}
\docType{methods}
\alias{hasTaxasets}
\alias{hasTaxasets-methods}
\alias{hasTaxasets,ANY-method}
\alias{hasTaxasets,list-method}
\alias{hasTaxasets,phylo4d-method}
\title{ Check if there are taxasets}
\description{
Check if there are any taxasets
}
\section{Methods}{
\describe{

\item{\code{signature(x = "ANY")}}{ Returns false, for any other object }

\item{\code{signature(x = "list")}}{ Checks the first object in the list for taxasets }

\item{\code{signature(x = "phylo4d")}}{ Check for taxasets }
}}
\keyword{methods}
