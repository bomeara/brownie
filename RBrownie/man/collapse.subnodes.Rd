\name{collapse.subnodes}
\alias{collapse.subnodes}
\title{
Collapse subnodes into singletons
}
\description{
Method which converts any existing subnodes into their phylobase equilvalent, singleton nodes (nodes with an in-degree and out-degree equal to one).
}
\usage{
collapse.subnodes(x, rm.ex.data = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of type phylo4d_ext which contains subnodse
}
  \item{rm.ex.data}{
Should data columns that only exist in subnodes.data be removed (TRUE is recommended).
}
}
\details{
This functions allows a user to convert from \code{RBrownie:::phylo4d_ext} objects to \code{phylobase:::phylo4d} objects without loosing any of the tree's topological information (or data).  This is in contrast to the \code{\link{collapse.singletons}} function which removes subnodes and returns a phylo4d tree without singletons or subnodes
}
\value{
An object of class \code{phylo4d} with singleton nodes where subnodes used to be
}
\author{
J. Conrad Stack
}

\seealso{
\code{\link{collapse.singletons}},
\code{\link{phyext}},
\code{\link{expand.singles}}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

