\name{checkMono}
\alias{checkMono}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check for monophyly
}
\description{
Check if the taxa in a specified taxaset form a monophyly in each of the trees in a list.
}
\usage{
checkMono(treeslist, taxaset, percent = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{treeslist}{
A list of phylogenetic trees (list of brownie objects)
}
  \item{taxaset}{
The taxaset to use (as column index)
}
  \item{percent}{
Should the results be returned as the number of trees in the set which contain this monophyly (\code{percent=FALSE}) or the percentage of trees (\code{percent=TRUE})
}
}
\details{
This function exists so that a set of trees can be checked to see if how many of the trees within contain the monophyly specified by taxaset.  To see a list of the taxasets available to test call \code{taxasets()}.
}
\value{
The number of trees (as percent or number) which contain the monophyly.
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
\code{\link{taxa.mono}},
\code{\link{taxa.charvect}},
\code{\link{which.mono}}
}
\examples{
data(parrot)
numbcontaining = checkMono(parrot,1,percent=FALSE)
numb = length(parrot)
numbcontaining / numb
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

