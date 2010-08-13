\name{getSubNodePosish}
\alias{getSubNodePosish}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Return the subnode position
}
\description{
This function returns the position of subnodes as a two-column matrix.  The two-columns are used for now because the position will eventually be represented as an interval to reflect the uncertainty of node placement.
}
\usage{
getSubNodePosish(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A phylo4d_ext
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A two-column matrix with the positions 
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
data(parrot)
getSubNodePosish(parrot[[1]])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

