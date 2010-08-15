\name{run.analysis}
\alias{run.analysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Execute a brownie-formatted nexus file
}
\description{
This is mainly a testing function.  It attempts to execute \code{filename} using the brownie core.  This function should be used instead of execute.brownie (probably)
}
\usage{
run.analysis(filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
%%     ~~Describe \code{filename} here~~
}
}
\details{
This function can be used to directly execute a brownie file.  It returns very raw text from the brownie core, but would be good for debugging or if you wanted to execute your own custom formatted file to test out brownie core functions that are not yet built in to RBrownie.
}
\value{
\item{treesout}{Newick text strings which the brownie core generated.  'discrete' when used with the option reconstruct=yes will generate these trees.  Often they are SIMMAP-formatted (v1.0)}
\item{textout}{Raw text of the analysis output.  Try using scan.textout and read.analysis.output to get this information into a more tractable data structure}


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
\code{\link{scan.textout}},
\code{\link{read.analysis.output}},
\code{\link{execute.brownie}}
}
\examples{
data(parrot)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

