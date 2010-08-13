\name{areTaxaPara}
\alias{areTaxaPara}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Do taxa form a paraphyletic clade
}
\description{
This function determines whether or not a set of taxa form a paraphyletic clade.
}
\usage{
areTaxaPara(x, taxind, with.respect.to = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A brownie object
 }
  \item{taxind}{
The column index of the taxaset to be checked
}
  \item{with.respect.to}{
The column index of a different taxaset to check against taxind: Is this taxaset paraphyletic with respect to the taxaset indicated by taxind. Experimental.
}
}
\details{
A group of taxa is said to be paraphyletic if the group contains its last common ancestor but does not contain all the descendants of that ancestor.  The the context of RBrownie, a paraphyletic clade is one where all the taxa not in the taxaset indicated by taxind themselves form a monophyletic clade.  See function definition for more details.
}
\value{
True if the taxa form a paraphyly, false otherwise.
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
taxasets(parrot)
areTaxaPara(parrot[[1]],1) # FALSE
areTaxaPara(parrot[[1]],2) # TRUE
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
