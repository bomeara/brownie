\name{addChoose}
\alias{addChoose}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add choose command to a brownie object (or list of brownie objects)
}
\description{
This function adds a character string to the the \code{commands} slot of a brownie object (or list of brownie objects) which instructs the brownie program to use a certain tree or continuous or discrete character in later analysis in the \code{commands} slot.
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
Index (integer) of the continuous character to use.
}
  \item{tree}{
Index (integer) of the tree to use.
}
  \item{discchar}{
Index (integer) of the discrete character to use.
}
}
\details{
The sQuote{choose} command instructs the brownie core to use a certain parameter for further analysis.   The \code{commands} slot is read sequentially by the Brownie program, so any sQuote{choose} commands should be placed above the analyses they are meant to effect.  Also note that this command is \emph{not} additive, meaning that consecutive sQuote{choose} commands using the same argument will not use both arguments, just the later.

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
#ptmp = addChoose(ptmp,discchar=2) # this would change the discchar from 1 to 2
commands(ptmp)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

