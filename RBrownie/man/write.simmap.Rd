\name{write.simmap}
\alias{write.simmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Write phylo4d_ext tree to newick string
}
\description{
This function write a phylo4d_ext object to a SIMMAP-formatted (version 1.0 or 1.1) text string.  This function is mainly used in conjunction with \code{write.nexus.simmap}.  It is an extension of the \code{write.tree} function in the ape package.
}
\usage{
write.simmap(x, usestate = NULL, file = "", vers = 1.1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A phylo4d_ext object
}
  \item{usestate}{
Which data column should be written in the internal simmap blocks.  If NULL, then either the first data column (vers==1.1 || vers==1.0) or all data columns (vers==1.5) are used.
}
  \item{file}{
The file to write to.  If not set, then the tree is printed to the screen
}
  \item{vers}{
Which version to use.  1.0 uses the ';' as the internal delimiter which 1.1 uses ':'
}
  \item{\dots}{
These are passed to write.tree (which is only called if there are no subnodes)
}
}
\details{

}
\value{
A character string (if \code{file} is not set)
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
\code{\link{newlabels.v15}},
\code{\link{newlabels.v1x}},
\code{\link{write.nexus.simmap}}
}
\examples{
data(parrot)
write.simmap(parrot[[1]]) # write to screen
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

