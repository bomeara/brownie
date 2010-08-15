\name{writeBrownie-methods}
\docType{methods}
\alias{writeBrownie-methods}
\alias{writeBrownie,brownie-method}
\alias{writeBrownie,list-method}
\title{ Write brownie object to file}
\description{
 This function writes out brownie objects (or lists of such) to a nexus file, which is slightly tailored so that it is compatible with the brownie core or the brownie command-line program.  It is used internally by RBrownie to write brownie objects to temporary file, which are then executed by the brownie core. 
}
\section{Methods}{
\describe{

\item{\code{signature(x = "brownie")}}{ This is single object is passed as a list to the 'list' definition of this function }

\item{\code{signature(x = "list")}}{  
Brownie object are mainly translated as standard NEXUS files with a few important exceptions.  First, the TREES block can contain newick tree strings which are in SIMMAP (v1.0 or v1.1) format.  The discrete and continuous data in the \code{data} slot are separated into two different blocks:  CHARACTERS and CHARACTERS2.  Any taxasets that exists are written to an ASSUMPTIONS block and anything in the \code{commands} slot is written within a BROWNIE block.
	
This function is still relatively untested and if any errors are returned please submit then to the author, who would greatly appreciate it.  
}
}}
\keyword{methods}
\keyword{ ~~ other possible keyword(s)}
