\name{checkPara}
\alias{checkPara}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check for paraphyly
}
\description{
Check if the taxa in a specified taxaset form a paraphyly in each of the trees in a list.
}
\usage{
checkPara(treeslist, taxaset, percent = T)
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
This function exists so that a set of trees can be checked to see if how many of the trees within contain the paraphyly specified by taxaset.  To see a list of the taxasets available to test call \code{taxasets()}.
}
\value{
The number of trees (as percent or number) which contain the paraphyly
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
numbcontaining = checkPara(parrot,2,percent=FALSE)  # check is the 2nd taxaset forms a paraphyly
numb = length(parrot)
numbcontaining / numb
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

