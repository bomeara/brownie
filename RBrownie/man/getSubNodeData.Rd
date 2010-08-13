\name{getSubNodeData}
\alias{getSubNodeData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Return subnode data
}
\description{
Return the subnode data.frame
}
\usage{
getSubNodeData(x, colname)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A phylo4d_ext object
}
  \item{colname}{
The column of the subnode data.frame to return
}
}
\details{

}
\value{
The data.frame or column
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
getSubNodeData(parrot[[1]])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

