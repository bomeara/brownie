\name{brownie-class}
\Rdversion{1.1}
\docType{class}
\alias{brownie-class}
\alias{brownie}
\alias{commands}
\alias{commands<-}
\alias{commands,brownie-method}
\alias{commands,list-method}
\alias{commands<-,brownie-method}
\alias{commands<-,list-method}
\alias{datatypes}
\alias{datatypes<-}
\alias{datatypes,brownie-method}
\alias{datatypes,list-method}
\alias{datatypes<-,brownie-method}
\alias{datatypes<-,list-method}
\alias{removeCommands}
\alias{removeCommands,list,ANY-method}
\alias{removeCommands,brownie,numeric-method}
\alias{removeCommands,brownie,character-method}
\alias{rmdata,brownie,numeric-method}
\alias{taxasets}
\alias{taxasets,phylo4d-method}
\alias{taxasets,brownie-method}
\alias{taxasets,list-method}
\alias{taxasets<-}
\alias{taxasets<-,brownie-method}
\alias{taxasets<-,list-method}
\alias{taxasets<-,phylo4d-method}
\alias{writeBrownie}
\alias{writeBrownie,brownie-method}
\alias{writeBrownie,list-method}

\title{Class "brownie" }
\description{ Brownie class }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("brownie", ...)} or by calls to brownie()

The 'brownie' S4 object builds off of the phylo4d_ext class, itself chiefly designed as an intermediary between phylo4d and brownie classes.  The brownie class should contain all the information necessary to run various brownie analyses.  

For that reason it supports subnodes (through \code{phylo4d_ext} objects) which are subnodes placed between nodes on a phylogenetic tree (and currently written out in SIMMAP v1.0 format).  It supports three different datatypes: taxasets (binary columns in the data slot which indicate membership in a taxaset with a value of 1), discrete data (data which you would want to represent as factors), and continuous data (numeric data).  It supports "commands", which are the actual text strings that are fed into the brownie core instructing it to carry out certain actions.  

RBrownie does not use a proprietary class to house multiple phylogenetic trees (which will be commonly needed).  Instead multiple brownie objects are put into a list and most of the functions below are also defined for the 'list' class.  Each brownie entry in the list contains the same information differing only in the which phylogenetic tree it contains.  Adding a taxaset or clearing the commands of a list of brownie object will therefore change every object in the list.  The hope is to fix this data redundancy in the future.

}
\section{Slots}{
  \describe{
    \item{\code{commands}:}{ Character string brownie commands }
    \item{\code{datatypes}:}{ Datatype values describing the data slot columns }
    \item{\code{subnode.id}:}{ Subnodes identifier.  Not really used }
    \item{\code{subnode.data}:}{ Data.frame which has the same structure as the data slot}
    \item{\code{subnode.branch}:}{ Which branch a subnode exists on }
    \item{\code{subnode.pos}:}{Where on the branch the subnode is positioned (relative to the descendant node as a fraction)}
    \item{\code{weight}:}{ The tree weight }
    \item{\code{data}:}{ Inherited from phylo4d }
    \item{\code{metadata}:}{ Inherited from phylo4d }
    \item{\code{edge}:}{ Inherited from phylo4d }
    \item{\code{edge.length}:}{ Inherited from phylo4d }
    \item{\code{label}:}{ Inherited from phylo4d }
    \item{\code{edge.label}:}{ Inherited from phylo4d }
    \item{\code{order}:}{ Inherited from phylo4d }
    \item{\code{annote}:}{ Inherited from phylo4d }
  }
}
\section{Extends}{
Class \code{"\linkS4class{phylo4d_ext}"}, directly.
Class \code{"\linkS4class{phylo4d}"}, by class "phylo4d_ext", distance 2.
Class \code{"\linkS4class{phylo4}"}, by class "phylo4d_ext", distance 3.
}
\section{Methods}{
  \describe{
    \item{clearCommands}{\code{signature(x = "brownie")}: 
    	Iterate through all brownie objects in this list and clear commands from each. 
    }
    \item{commands}{\code{signature(x = "brownie")}: ... }
    \item{commands<-}{\code{signature(x = "brownie")}:
    	options: add=FALSE,index=NULL,replace=FALSE
    	This function adds or appends commands to a brownie object or list of brownie object.
    }
    \item{datatypes}{\code{signature(x = "brownie")}: Return datatypes}
    \item{datatypes<-}{\code{signature(x = "brownie")}: Assign datatypes }
    \item{removeCommands}{\code{signature(x = "brownie", index = "numeric")}: remove certain brownie commands based on the row index }
    \item{removeCommands}{\code{signature(x = "brownie", index = "character")}: remove certain brownie commands which start with 'index'. }
    \item{rmdata}{\code{signature(x = "brownie", index = "numeric")}: remove the data given the index }
    \item{taxasets}{\code{signature(x = "brownie")}: display taxasets }
    \item{taxasets<-}{\code{signature(x = "brownie")}: This function adds a taxaset.  Taxaset names are appended with the string TAXASET_.  This is one way brownie can tell if a certain column in the data slot is a taxaset. }
    \item{writeBrownie}{\code{signature(x = "brownie")}: 
 This function writes out brownie objects (or lists of such) to a nexus file, which is slightly tailored so that it is compatible with the brownie core or the brownie command-line program.  It is used internally by RBrownie to write brownie objects to temporary file, which are then executed by the brownie core. 

Brownie object are mainly translated as standard NEXUS files with a few important exceptions.  First, the TREES block can contain newick tree strings which are in SIMMAP (v1.0 or v1.1) format.  The discrete and continuous data in the \code{data} slot are separated into two different blocks:  CHARACTERS and CHARACTERS2.  Any taxasets that exists are written to an ASSUMPTIONS block and anything in the \code{commands} slot is written within a BROWNIE block.  This function is still relatively untested and if any errors are returned please submit then to the author, who would greatly appreciate it. 
 
    }
	 }
}
\references{}
\author{ J. Conrad Stack}
\note{}

\seealso{

}
\examples{
showClass("brownie")


# remove commands:
data(geospiza_ext)
commands(geospiza_ext)
commands(removeCommands(geospiza_ext,"discrete")) # remove all 'discrete' commands
commands(clearCommands(geospiza_ext)) # remove all commands
}
\keyword{classes}
