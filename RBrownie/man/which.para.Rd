\name{which.para}
\alias{which.para}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find trees in a list which contain a paraphyly for a certain taxaset
}
\description{
This looks through a list of trees and returns the indices where the tree contains a paraphyly for the taxaset specified
}
\usage{
which.para(treeslist, taxaset)
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
A vector of indices of the list where the trees contained this paraphyly
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
\code{\link{which.mono}}
}
\examples{
data(parrot)
which.para(parrot,2) # returns empty vector
which.para(parrot,1) # all trees contain this paraphyly

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

