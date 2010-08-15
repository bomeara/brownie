\name{which.mono}
\alias{which.mono}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find trees in a list which contain a monophyly for a certain taxaset
}
\description{
This looks through a list of trees and returns the indices where the tree contains a monophyly for the taxaset specified
}
\usage{
which.mono(treeslist, taxaset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{treeslist}{
A list of brownie objects
}
  \item{taxaset}{
The taxaset index to check
}
}
\details{
This command is similar to functions likt \code{which} and \code{which.min}.
}
\value{
A vector of indices of the list where the trees contained this monophyly
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
\code{\link{which.para}}
}
\examples{
data(parrot)
which.mono(parrot,2) # returns empty vector
which.mono(parrot,1) # all trees contain this monophyly
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

