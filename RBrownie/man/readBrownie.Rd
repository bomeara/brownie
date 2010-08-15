\name{readBrownie}
\alias{readBrownie}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read in a brownie formatted nexus file
}
\description{
This function reads a brownie-formatted nexus file into a brownie object.  The nexus blocks it reads are: TAXA, TREES, CHARACTERS, CHARACTERS2, ASSUMPTIONS, and BROWNIE. 
}
\usage{
readBrownie(fname)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fname}{
%%     ~~Describe \code{fname} here~~
}
}
\details{
The CHARACTERS and CHARACTERS2 blocks are both read into the data slot.  The ASSUMPTIONS block is scanned for the string 'taxset name = <list of taxa>', it then sets up another column in the data slot for this taxaset.  All the text from the BROWNIE block is put into the commands slot - it is not checked for accuracy or valid commands at this point.

If more than one tree is found in the TREES block, then multiple brownie objects are created (one for each tree) and returned together as a list.  All of the data is duplicated so that each part of the list is a self-contained brownie object.
}
\value{
A brownie object or list of brownie objects.
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
\dontrun{readBrownie("parrot.nex")}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

