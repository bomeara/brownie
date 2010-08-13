\name{check.empty}
\alias{check.empty}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check character vectors
}
\description{
Check for empty strings in a character vector
}
\usage{
check.empty(bb)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bb}{
Character vector
}
}
\details{
Mainly an internal function.
}
\value{
Vector of logicals indicating if a line is an empty string
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
junk = c("asdfasd","","asdfasdf")
check.empty(junk)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

