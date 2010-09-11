\name{addSubNode}
\alias{addSubNode}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add subnode to \code{phylo4d_ext} object
}
\description{
Add a subnode to \code{phylo4d_ext} object.  
}
\usage{
addSubNode(x, anc, dec, position, dataf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
is a phylo4d_ext object 
}
  \item{anc}{
Node id of the ancestor on the branch on which the subnode is to be added
}
  \item{dec}{
Node id of the descendant on the branch on which the subnode is to be added
}
  \item{position}{
The position on the branch where the subnode is to be placed.  Position is currently a 2 column matrix - this is because it will eventually be used to store confidence intervals.  For now it just the proportion of the branch starting from the ancestor.
}
  \item{dataf}{
a data.frame object which mimics the structure of the the data slot.  It contains the subnode data (which is why it should look similar to data).
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
\code{\link{phyext}}
}
\examples{


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
