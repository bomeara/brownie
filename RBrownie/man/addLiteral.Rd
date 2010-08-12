\name{addLiteral}
\alias{addLiteral}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add command string to brownie object
}
\description{
Add command string to brownie object manually
}
\usage{
addLiteral(obj, cmdstr = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
A brownie object (or list of brownie objects)
}
  \item{cmdstr}{
A character string which is added to commands.
}
}
\details{
Caution should be used when using this function.  It adds cmdstr to the commands slot of a brownie object.  It adds the character string as is, so this command really should only be used by those people who are familiar with brownie command-line commands.
}
\value{
Returns the original object (obj) with ratetest command added.
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
ptmp = clearCommands(parrot)
# this will probably crash the brownie core if actually executed
ptmp = addLiteral(ptmp,"Anything goes!") 
ptmp = addLiteral(ptmp,"...but not really") 
commands(ptmp)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
