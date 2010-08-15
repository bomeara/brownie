\name{taxa.charvect}
\alias{taxa.charvect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Return a character vector of taxa
}
\description{
Convert a taxaset index (column in the data frame) to a character vector containing the names of the taxa. 
}
\usage{
taxa.charvect(x, taxind, append.internal = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A brownie object which contains taxa sets
}
  \item{taxind}{
The column index of the taxaset to transcribe, or a character string specifying the name of the taxaset
}
  \item{append.internal}{
For internal use only. 
}
}
\details{
This functions finds the taxaset column in the \code{data} slot of a \code{x}.  It then returns the tipLabels of the taxa where there is a 1 in the data column.  (Brownie uses binary tip data to define a taxaset: 0 - not a member of the taxaset; 1 = member).
}
\value{
A character vector containing the names of the taxa in this set
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
taxa.charvect(parrot[[2]],'intrajoint')
taxa.charvect(parrot[[1]],2) # should be the same thing

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

