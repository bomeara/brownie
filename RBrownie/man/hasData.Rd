\name{hasData-methods}
\docType{methods}
\alias{hasData}
\alias{hasData-methods}
\alias{hasData,list-method}
\alias{hasData,phylo-method}
\alias{hasData,phylo4-method}
\alias{hasData,phylo4d-method}
\title{ Methods for Function hasData in Package `RBrownie' }
\description{
This method checks to see if tree objects have data
}
\section{Methods}{
\describe{

\item{\code{signature(x = "list")}}{
Calls hasData again on each element of the list
}

\item{\code{signature(x = "phylo")}}{
Returns FALSE
}

\item{\code{signature(x = "phylo4")}}{
Returns FALSE
}

\item{\code{signature(x = "phylo4d")}}{
Check if \code{ncol(tdata(x))} is equal to zero or not
}
}}
\keyword{methods}

