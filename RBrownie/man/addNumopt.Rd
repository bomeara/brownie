\name{addNumopt}
\alias{addNumopt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add numerical optimization command to brownie object
}
\description{
This function adds a character string to the the commands slot of a brownie object (or list of brownie objects) which instructs the brownie program to set certain variables when performing numerical optimization.
}
\usage{
addNumopt(obj, 
	iter = NULL, 
	toler = NULL, 
	randstart = NULL, 
	seed = NULL, 
	stepsize = NULL, 
	detail = FALSE, 
	redo = FALSE, 
	giveupfactor = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
A brownie object (or list of brownie objects)
}
  \item{iter}{
sets the maximum number of iterations of the Nelder-Mead simplex algorithm
}
  \item{toler}{
Sets the precision of the stopping criterion: what amount of change in the likelihood is considered small enough to count as zero change.
}
  \item{randstart}{
Sets the number of random starts to use.
}
  \item{seed}{
Sets the seed for the random variable
}
  \item{stepsize}{
sets the NM step size.
}
  \item{detail}{
Specifies whether or not to have detailed output from numerical optimization
}
  \item{redo}{
Specifies whether to redo reps which stop due to iteration limits
}
  \item{giveupfactor}{
When \code{redo=T}, is used to tell the software when to stop restarting: when the ratio of unsuccessful to successful starts is > \code{giveupfactor}.
}
}
\details{
Sets options for numerical optimization.
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

## The function is currently defined as
data(parrot)
ptmp = clearCommands(parrot)
ptmp = addNumopt(ptmp,iter=1000,randstart=50)
commands(ptmp)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
