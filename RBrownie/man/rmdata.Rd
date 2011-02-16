\name{rmdata-methods}
\docType{methods}
\alias{rmdata}
\alias{rmdata-methods}
\alias{rmdata,brownie,character,missing-method}
\alias{rmdata,brownie,numeric,ANY-method}
\alias{rmdata,brownie,numeric,missing-method}
\alias{rmdata,brownie,numeric,numeric-method}
\alias{rmdata,list,ANY,missing-method}
\alias{rmdata,phylo4d,character,missing-method}
\alias{rmdata,phylo4d_ext,character,missing-method}
\alias{rmdata,phylo4d_ext,numeric,missing-method}
\alias{rmdata,phylo4d_ext,numeric,numeric-method}
\alias{rmdata,phylo4d,numeric,missing-method}
\title{ rmdata from a tree object}
\description{
This is the recommended way of removing data from the \code{data} and \code{subnode.data} slots of 'phylo4d_ext' and 'brownie' objects.  Although this method is defined for both 'brownie' and 'phylo4d_ext' classes, there is not much difference between them.
}
\section{Methods}{
\describe{

\item{\code{signature(x = "brownie", index = "character", subindex = "missing")}}{
Enter the column name of the data column to be removed (see \code{colnames(tdata(x))} for valid options). \code{subindex} is set to the value of \code{index}. This function will attempt to remove that column from both regular and subnodes.
}

\item{\code{signature(x = "brownie", index = "numeric", subindex = "ANY")}}{
Enter the column index of the data column to be removed.  When \code{subindex} is ANY, it is ignored (except in the cases below) and set to the value of \code{index}
}

\item{\code{signature(x = "brownie", index = "numeric", subindex = "missing")}}{
Enter the column index of the data column to be removed (\code{seq(ncol(tdata(x)))} are valid options). \code{subindex} is set to the value of \code{index}.  This function will attempt to remove that column from both regular and subnodes, regardless of whether their column names match up.
}

\item{\code{signature(x = "brownie", index = "numeric", subindex = "numeric")}}{
Used internally: remove the data column (\code{index}) and the subnode data column (\code{subindex}) from their respective slots
}

\item{\code{signature(x = "list", index = "ANY", subindex = "missing")}}{
Applys \code{rmdata} to each object in a list
}

\item{\code{signature(x = "phylo4d", index = "character", subindex = "missing")}}{
Enter the column name of the data column to be removed.  No need to specify a subindex
}

\item{\code{signature(x = "phylo4d_ext", index = "character", subindex = "missing")}}{
Enter the column name of the data column to be removed (see \code{colnames(tdata(x))} for valid options). \code{subindex} is set to the value of \code{index}. This function will attempt to remove that column from both regular and subnodes.
}

\item{\code{signature(x = "phylo4d_ext", index = "numeric", subindex = "missing")}}{
Enter the column index of the data column to be removed (\code{seq(ncol(tdata(x)))} are valid options). \code{subindex} is set to the value of \code{index}. This function will attempt to remove that column from both regular and subnodes, regardless of whether their column names match up.
}

\item{\code{signature(x = "phylo4d_ext", index = "numeric", subindex = "numeric")}}{
Used internally: remove the data column (\code{index}) and the subnode data column (\code{subindex}) from their respective slots
}

\item{\code{signature(x = "phylo4d", index = "numeric", subindex = "missing")}}{
Enter the column index of the data column to be removed.  No need to specify a subindex
}
}}
\keyword{methods}

