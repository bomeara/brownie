\name{brownie-methods}
\docType{methods}
\alias{brownie-methods}
\alias{brownie,list-method}
\alias{brownie,phylo-method}
\alias{brownie,phylo4-method}
\alias{brownie,phylo4d-method}
\title{ Methods for generic brownie constructor function }
\description{
 These functions all create brownie objects out of objects that the brownie class extends (or is related to).
}
\section{Methods}{
\describe{

\item{\code{signature(x = "list")}}{ Convert a list of any of the below class types into a list of brownie objects }

\item{\code{signature(x = "phylo")}}{ Create brownie object from an object of ape's phylo class.  The object is first coersed into a phylo4 object (using phylobase methods) and then passed to the phylo4 brownie-constructor. }

\item{\code{signature(x = "phylo4")}}{ Create brownie object from an object of phylobase's phylo4 S4 class.  The object is converted into a phylo4d object (using phylobase methods) with no data and then passed to the phylo4d brownie-constructor.  }

\item{\code{signature(x = "phylo4d")}}{ Calls \code{new} function to create a brownie object.  If there is data in the \code{data} slot of the phylo4d object, then brownie attempts to guess the correct datatype. 

There is an argument that can be passed to this function called \code{known.types} where you can specify which datatypes you want to use for each column in the \code{data} slot's data.frame.   

See the RBrownie:::.guess.datatype function for how brownie guesses datatypes and the \code{\link{brownie.datatypes}} function to enumerate the valid datatypes.
}
}}
\keyword{methods}
\keyword{ ~~ other possible keyword(s)}
