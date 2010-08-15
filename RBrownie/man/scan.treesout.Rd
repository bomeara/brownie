\name{scan.treesout}
\alias{scan.treesout}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Scan raw brownie tree text
}
\description{
This is mainly an internal function.  It scans raw text returned from brownie (in the treesout)
}
\usage{
scan.treesout(output)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{output}{
The brownie character string
}
}
\details{
This scans the treesout section of brownie output for SIMMAP formatted newick strings. It then converts them to phylo4d_ext objects
}
\value{
A list of phylo4d_ext objects
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
data(geospiza_ext)
rawtext = run.asis(geospiza_ext)
scan.treesout(rawtext$treesout)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

