\name{read.continuous.output}
\alias{read.continuous.output}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read brownie output from Noncensored analysis
}
\description{
This function processes output returned from the brownie core before passing it on to read.analysis.output.
}
\usage{
read.continuous.output(filename, txt = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
A file with brownie output
}
  \item{txt}{
Text containing brownie output
}
  \item{\dots}{
Not currently used
}
}
\details{
This function removes WARNINGS returned by the brownie core output before passing the rest of it on to read.analysis.output.  This is also somewhat of an internal function
}
\value{
A data.frame with the discrete analysis results
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
\code{\link{read.analysis.output}}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

