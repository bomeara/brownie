\name{scan.textout}
\alias{scan.textout}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Scan raw brownie analysis text
}
\description{
This is mainly an internal function.  It scans raw text returned from brownie (in the textout)
}

\usage{
scan.textout(output, strip.white = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{output}{
The brownie character string
}
  \item{strip.white}{
Should the white space should be striped from either end of the character string?
}
}
\details{
This function also removes any empty string lines
}
\value{
A character vector
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
\dontrun{
data(parrot)
junkrun=run.asis(parrot)
cuttxt = scan.textout(junkrun$textout[1]) # post-process the text
}
}

\keyword{ methods }

