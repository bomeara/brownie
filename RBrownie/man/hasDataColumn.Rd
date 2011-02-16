\name{hasDataColumn-methods}
\docType{methods}
\alias{hasDataColumn}
\alias{hasDataColumn-methods}
\alias{hasDataColumn,phylo4,ANY-method}
\alias{hasDataColumn,phylo4d,character-method}
\alias{hasDataColumn,phylo4d,numeric-method}
\alias{hasDataColumn,phylo,ANY-method}
\title{ Does a tree object contain a specific data column }
\description{
This function is used to query whether or not a class (inheriting from 'phylo4d') contains a specific data column defined either by its column index or it's column name
}
\section{Methods}{
\describe{

\item{\code{signature(x = "phylo4", index = "ANY")}}{
Return FALSE
}

\item{\code{signature(x = "phylo4d", index = "character")}}{
Checks whether or not the data slot contains a column with the name \code{index} 
}

\item{\code{signature(x = "phylo4d", index = "numeric")}}{
Checks whether or not \code{index} is less than or equal to \code{ncol(tdata(x))}
}

\item{\code{signature(x = "phylo", index = "ANY")}}{
Return FALSE
}
}}
\keyword{methods}

