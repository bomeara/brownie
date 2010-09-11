\name{addStartLog}
\alias{addStartLog}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add log start command to brownie object
}
\description{
This function adds a character string to the the commands slot of a brownie object (or list of brownie objects) which instructs the brownie program to start logging to a certain file.  (Be careful when using this.  Logging doesn't seem stable in the brownie core).
}
\usage{
addStartLog(obj, 
	file = NULL, 
	fileappend = FALSE, 
	filereplace = FALSE, 
	usetempfile = FALSE, 
	...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
A brownie object (or list of brownie objects)
}
  \item{file}{
File to log to
}
  \item{fileappend}{
Should File be appended to?
}
  \item{filereplace}{
Should File be replaced?
}
  \item{usetempfile}{
Should a tempfile be written to (instead of file).  Not currently used.
}
  \item{\dots}{
Not currently used.
}
}
\details{
Adds 'log start <options>;' 
}
\value{
Returns the original object (obj) with start command added.
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
ptmp = addStartLog(ptmp,file="BrownieLog.txt",filereplace=TRUE)
commands(ptmp)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
