\name{removeTaxasets-methods}
\docType{methods}
\alias{removeTaxasets}
\alias{removeTaxasets-methods}
\alias{removeTaxasets,brownie,character-method}
\alias{removeTaxasets,brownie,numeric-method}
\title{ Remove taxasets from a brownie object }
\description{
This function attempts to remove taxaset data columns from a brownie object by matching a taxaset to a column and then removing that column with \code{rmdata}.
}
\section{Methods}{
\describe{

\item{\code{signature(x = "brownie", index = "character")}}{
Remove the taxaset named \code{index}
}

\item{\code{signature(x = "brownie", index = "numeric")}}{
Remove the \code{index}-th taxaset. 
}
}}
\keyword{methods}

