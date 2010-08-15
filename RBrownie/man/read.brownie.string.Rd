\name{read.brownie.string}
\alias{read.brownie.string}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read brownie command string into a list
}
\description{
Mainly an internal function.  It tokenizes a character string - might be used in the future to check whether or not a brownie command is valid.
}
\usage{
read.brownie.string(txtstr, operator = "=")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{txtstr}{
%%     ~~Describe \code{txtstr} here~~
}
  \item{operator}{
%%     ~~Describe \code{operator} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  \item{command}{The main brownie command}
%%  \item{options}{The options used with the main command}
%% ...
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
read.brownie.string(commands(parrot)[3])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

