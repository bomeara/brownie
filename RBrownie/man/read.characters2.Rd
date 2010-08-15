\name{read.characters2}
\alias{read.characters2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read CHARACTERS2 block
}
\description{
This function reads raw text in from the nexus file finput.  It then writes the TAXA and \code{blockname} blocks to a temporary nexus file using CHARACTERS instead of \code{blockname}.  It then re-reads the temporary nexus file using \code{readNexus(tmpfile,type="data")}.  
}
\usage{
read.characters2(finput, blockname = "characters2")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{finput}{
%%     ~~Describe \code{finput} here~~
}
  \item{blockname}{
%%     ~~Describe \code{blockname} here~~
}
}
\details{
This function was created mainly to handle reading CHARACTERS2 nexus blocks.  These blocks are used by brownie along with CHARACTERS to separate discrete (factor) data from continuous (numeric) data.  It is also used because SIMMAP-formatted (version < 1.5) nexus files tend to cause problems with the nexus class library (which phylobase uses)
}
\value{
A data.frame containing the data from \code{blockname}.
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
\code{\link{readNexus}}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

