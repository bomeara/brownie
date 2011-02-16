\name{read.nexus.simmap}
\alias{read.nexus.simmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read a nexus file with SIMMAP-formatted trees (v1.0)
}
\description{
This function reads the trees from a nexus file where the trees are SIMMAP-formatted (v1.0).  It does not read any data blocks, just the trees.  If the file contains multiple trees, this function returns them as a list.
}
\usage{
read.nexus.simmap(finput = "", text = NULL, vers = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{finput}{
A nexus file
}
  \item{text}{
Alternatively, read this text.
}
  \item{vers}{
Which version of simmap to try to read.
}
  \item{\dots}{
All other arguments are passed to internal \code{read.simmap} function.
}
}
\details{
This function reads nexus blocks using a combination of internal functions and a few phylobase functions.  If simmap newick trees are found then the function \code{read.simmap} is used otherwise a modified \code{read.tree} is found.
}
\value{
A list containing either phylo4d objects or phylo4d_ext objects depending on whether or not SIMMAP characters are found.
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

