\name{collapse.singletons}
\alias{collapse.singletons}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Collapse singleton nodes 
}
\description{
This function turns any singleton nodes (internal non-root nodes with degree = 2) into subnodes.  It uses the function \code{addSubNodes} to attach the nodes to a phylo4d_ext object.  These subnodes are positioned on the branch relative to the descendant node.  
}
\usage{
collapse.singletons(phy)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{
phylo4 tree with singletons that will potentially be collapsed.
}
}
\details{

}
\value{
The tree with the singletons collapsed.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
J. Conrad Stack
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{hasSubNodes}}
}
\examples{

# example text:
phytext = "((Taxon1:{A,0.1; C,0.1}, Taxon2:{T,0.1; C,0.1}):{C,0.5}, Taxon3:{C,0.4} );"
# read text as a simmap tree
phy = read.simmap(text=phytext,vers=1.0) 
hasSingle(phy) # TRUE
phy = collapse.singletons(phy)
hasSingle(phy) # FALSE

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

