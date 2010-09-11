\name{addChoose}
\alias{addChoose}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add choose command to a brownie object (or list of brownie objects)
}
\description{
This function adds a character string to the the \code{commands} slot of a brownie object (or list of brownie objects) which instructs the brownie program to "choose" or select a certain tree, continuous or discrete character for further analysis.
}
\usage{
addChoose(obj, char = NULL, tree = NULL, discchar = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
A brownie object (or list of brownie objects)  
}
  \item{char}{
Choose a certain continuous character for analysis.
}
  \item{tree}{
Choose a certain tree character for analysis.
}
  \item{discchar}{
Choose a certain discrete character for analysis.
}
}
\details{
The 'choose' command instructs the brownie core to use a certain parameter for further analysis.  Only use one of the options.

From brownie core:
Selects one tree and/or character to use for subsequent analyses. Continuous characters are chosen by default ('char=' or 'c='), discrete characters must be chosen using 'discretechar=' or 'd='.

The remaining trees and characters are still stored in memory and can be chosen later.

}
\value{
Returns the original object (\code{obj}) with ratetest command added to its \code{commands} slot.
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

# These will rarely be chained like this, but for example:
ptmp = addChoose(ptmp,char=5)
ptmp = addChoose(ptmp,tree=2)
ptmp = addChoose(ptmp,discchar=1)
commands(ptmp)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

