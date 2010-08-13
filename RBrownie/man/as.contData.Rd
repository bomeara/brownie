\name{as.contData}
\alias{as.contData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert to continuous data
}
\description{
Essentially this function will convert factor data to numeric data, which is equivalent to the distinction that brownie makes between continuous and discrete data.
}
\usage{
as.contData(dat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
  A vector of undefined type
}
}
\details{
Brownie uses the factor and numeric types to determine whether or not a data column in the data slot of a brownie object is 'discrete' or 'continuous'.  The distinction between the two datatypes is important as the two types are written to two different nexus blocks when brownie writes to file (CHARACTERS and CHARACTERS2).

This is mainly an internal function, but could be useful as a way to make sure your data is compatible with brownie datatypes (which might change in the future, but this function call will not)
}
\value{
Returns dat converted to a numeric vector
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
x = sample(1:3,10,replace=TRUE)
class(x)
x = as.discData(x)
class(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
