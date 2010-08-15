\name{showSubNodes}
\alias{showSubNodes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Display subnodes
}
\description{
Display a phylo4d_ext object.  This mainly extends the show function for phylo4 object in that it also displays any subnodes if they exist.
}
\usage{
showSubNodes(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A phylo4d_ext object
}
}
\details{
Attempts to print out the branches which contain subnodes and where on the branches the subnodes are.  This is still a bit experimental.
}
\value{
Prints out the phylo4d_ext object.
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
showSubNodes(parrot[[1]])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

