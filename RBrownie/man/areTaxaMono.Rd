\name{areTaxaMono}
\alias{areTaxaMono}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Do taxa form a monophyletic clade
}
\description{
Check whether or not a set of taxa (existing withing taxasets()) together form a monophyletic clade.
}
\usage{
areTaxaMono(x, taxindex)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A brownie object
}
  \item{taxindex}{
The column index of the taxaset to be checked
}
}
\details{
Use the command \code{taxasets(x)} to see which taxasets are available for checking.  If the taxindex entered is greater than the number of taxasets then this command simply returns false and a warning.  The \code{taxasets(x)} command can also be assigned to in order to create new taxasets.
}
\value{
True if the taxa form a monophyly, false otherwise.
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
\code{\link{taxasets}}
}
\examples{
data(parrot)
taxasets(parrot)
areTaxaMono(parrot[[1]],1) # TRUE
areTaxaMono(parrot[[1]],2) # FALSE
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

