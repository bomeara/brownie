\name{get.left.right}
\alias{get.left.right}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get right and left tokens
}
\description{
Get character tokens to the left and right of char
}
\usage{
get.left.right(str, char)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{str}{
A character string to be tokenized.
}
  \item{char}{
Split the string on this character
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The two strings to the left and right of the split char.
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
get.left.right("&W=-12.52","=")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

