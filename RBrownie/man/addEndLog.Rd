\name{addEndLog}
\alias{addEndLog}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add log stop command to brownie object
}
\description{
This function adds a character string to the the \code{commands} slot of a brownie object (or list of brownie objects) which instructs the brownie program to stop logging to a certain file.  (Be careful when using this.  Logging doesn't seem stable in the brownie core).

}
\usage{
addEndLog(obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
A brownie object (or list of brownie objects)
}
}
\details{
Adds 'log stop;' to the \code{commands} slot of a brownie object
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Returns the original object (obj) with ratetest command added.
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
