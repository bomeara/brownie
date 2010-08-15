\name{phylo4d_ext-class}
\Rdversion{1.1}
\docType{class}
\alias{phylo4d_ext-class}
\alias{hasSubNodes,phylo4d_ext-method}
\alias{hasWeight,phylo4d_ext-method}
\alias{plot,phylo4d_ext,missing-method}
\alias{rmdata,phylo4d_ext,numeric-method}
\alias{snbranch,phylo4d_ext-method}
\alias{sndata,phylo4d_ext-method}
\alias{sndata<-,phylo4d_ext-method}
\alias{snid,phylo4d_ext-method}
\alias{snposition,phylo4d_ext-method}
\alias{weight,phylo4d_ext-method}
\alias{weight<-,phylo4d_ext-method}

\title{Class "phylo4d_ext" ~~~ }
\description{ This class extends phylo4d in adding support for subnodes (not singletons) and tree weights.  }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("phylo4d_ext", ...)}.
	 ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{subnode.id}:}{Object of class \code{"integer"}. The subnode identifier.  Not currently being used. }
    \item{\code{subnode.data}:}{Object of class \code{"data.frame"}.  The subnode data frame.  This data frame should extend the parent phylo4d class data frame - it extends the parent data frame in the same way the phylo4d_ext extends phylo4d. }
    \item{\code{subnode.branch}:}{Object of class \code{"matrix"}.  This is a two-column matrix.  The first column contains the ancestor node Id and the descendant node Id for the two nodes spanning a branch.  This copies data directly from the \code{edge} slot.  It uses the node ids instead of the edge row index because the edge slot is often reordered.  }
    \item{\code{subnode.pos}:}{Object of class \code{"matrix"}.  This is the position of the subnode on the branch as a fraction of the branch length.  It is the fraction of the branch length starting from the descendant node. 
    	
    	So for the hypothetical branch: ANC <--------subnode--> DEC, the subnode position might be 0.25.
    	
    	This is currently a two column matrix, but the columns should contain the same data.  Eventually, these two columns will represent a confidence interval that the subnode might exist within.
     }
    \item{\code{weight}:}{Object of class \code{"numeric"}.  This is the weight for this tree.  The weight is gleaned from nexus comments in a nexus file or set manually.
    
    	This is the only piece of data that is not replicated across every member of a phylo4d_ext list.  Each tree can have a unique weight.
    }
    \item{\code{data}:}{Object of class \code{"data.frame"}.  Inherited from phylo4d }
    \item{\code{metadata}:}{Object of class \code{"list"}.  Inherited from phylo4d }
    \item{\code{edge}:}{Object of class \code{"matrix"}.  Inherited from phylo4 }
    \item{\code{edge.length}:}{Object of class \code{"numeric"}.  Inherited from phylo4 }
    \item{\code{label}:}{Object of class \code{"character"}.  Inherited from phylo4 }
    \item{\code{edge.label}:}{Object of class \code{"character"}.  Inherited from phylo4 }
    \item{\code{order}:}{Object of class \code{"character"}.  Inherited from phylo4 }
    \item{\code{annote}:}{Object of class \code{"list"}.  Inherited from phylo4 }
  }
}
\section{Extends}{
Class \code{"\linkS4class{phylo4d}"}, directly.
Class \code{"\linkS4class{phylo4}"}, by class "phylo4d", distance 2.
}
\section{Methods}{
  \describe{
    \item{hasSubNodes}{\code{signature(x = "phylo4d_ext")}: Check whether or not a phylo4d_ext object or list of these objects contain subnodes.  Basically, checks the length of the subnode.id slot. }
    \item{hasWeight}{\code{signature(x = "phylo4d_ext")}:  Checks whether or not a phylo4d_ext object or list of these objects contain weights. }
    \item{plot}{\code{signature(x = "phylo4d_ext", y = "missing")}: Calls \code{phyextPlot} with no extra arguments }
    \item{rmdata}{\code{signature(x = "phylo4d_ext", index = "numeric")}: Remove data from a phylo4d_ext object.  If subnode.data exists, then the same data index is removed from this as well. }
    \item{snbranch}{\code{signature(x = "phylo4d_ext")}: Return the subnode.branch matrix }
    \item{sndata}{\code{signature(x = "phylo4d_ext")}: Return the subnode.data data.frame }
    \item{sndata<-}{\code{signature(x = "phylo4d_ext")}:  This adds column to the subnode data.frame.  Pass the datnames option which will be the column name(s) of the new data columns. }
    \item{snid}{\code{signature(x = "phylo4d_ext")}: Return the subnode.id vector }
    \item{snposition}{\code{signature(x = "phylo4d_ext")}: Return the subnode.pos matrix }
    \item{weight}{\code{signature(x = "phylo4d_ext")}: Return the weight for a phylo4d_ext object or a vector of weights for a list of phylo4d_ext objects. }
    \item{weight<-}{\code{signature(x = "phylo4d_ext")}: Replaces the weight vector. }
	 }
}
\references{ ~put references to the literature/web site here ~ }
\author{ J. Conrad Stack }
\note{ ~~further notes~~ }

\seealso{
	\code{\linkS4class{phylo4d}}
}
\examples{
showClass("phylo4d_ext")
data(parrot)
newobj = as(parrot[[1]],'phylo4d_ext')
}
\keyword{classes}
