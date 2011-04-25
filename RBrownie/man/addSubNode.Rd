\name{addSubNode}
\alias{addSubNode}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add subnode to phylo4d_ext object
}
\description{
Add a subnode to phylo4d_ext object.  
}
\usage{
addSubNode(x, anc, dec, position, dataf, pos.is.fraction = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a phylo4d_ext object
}
  \item{anc}{
Node id of the ancestor on the branch on which the subnode is to be added
}
  \item{dec}{
Node id of the descendant on the branch on which the subnode is to be added
}
  \item{position}{
The position on the branch where the subnode is to be placed.  It is stored as a faction of the branches' total length and is in relation to the ancestral node.  In other words, a smaller value will be closer to the ancestral node.
}
  \item{dataf}{
a data.frame object which should be structurally similar to the data frame of the phylo4d object (data slot).  It should contain data entries for the subnode.
}
  \item{pos.is.fraction}{
Specify whether \code{position} argument a fraction of the branch length or an absolute number.
}
}
\details{
This is the main function in RBrownie for adding subnodes manually.  It is used internally when reading in SIMMAP-formatted newick trees and can also be used to edit phylo4d_ext trees.
}
\value{
Returns \code{x} with a new subnode added.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{phyext}},
\code{\link{phyextPlot}}
}
\examples{


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
