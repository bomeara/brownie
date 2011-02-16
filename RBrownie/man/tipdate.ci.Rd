\name{tipdate.ci}
\alias{tipdate.ci}

\title{
Get coalescent intervals for a tip-dated tree
}
\description{
Similar to ape's coalescent intervals, except it works when the tree is not ultrametric.  This is being used for now as a helper function for treeHeight
}
\usage{
tipdate.ci(tr, show.plot = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tr}{
A tree in phylo format
}
  \item{show.plot}{
Should a plot be generated highlighting the coalescent events
}
}
\details{
This is mainly being used with treeHeights for now, but it is a good way to get the coalescent intervals for a tree.  Be aware though that the tree needs to be down-convert to \code{ape:::phylo} format first.
}
\value{
\item{lineages}{The number of extant lineages during this interval}
\item{interval.length}{The length of each interval}
\item{interval.count}{The total number of coalescent intervals}
\item{total.depth}{The height or depth of the tree}
\item{I}{Indicator for whether or not the event is a coalescent event (1) or a sampling event (0)}
}

\author{
J. Conrad Stack
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
tr = rtree(10)
ci = tipdate.ci(tr)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

