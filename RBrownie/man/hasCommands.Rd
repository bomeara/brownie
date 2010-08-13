\name{hasCommands}
\alias{hasCommands}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check for commands
}
\description{
Check if a brownie object or list of brownie objects contains any commands
}
\usage{
hasCommands(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A brownie object or list of objects
}
}
\details{
This function only checks if any commands exist
}
\value{
True if there are commands
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
hasCommands(parrot)
hasCommands(parrot[[1]])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

