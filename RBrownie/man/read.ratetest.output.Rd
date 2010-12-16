\name{read.ratetest.output}
\alias{read.ratetest.output}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read brownie output from Censored ratetest analysis
}
\description{
This function processes output returned from the brownie core before passing it on to \code{read.analysis.output}.
}
\usage{
read.ratetest.output(filename, txt = NULL, ...)
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
This function removes WARNINGS returned by the brownie core output before passing the rest of it on to read.analysis.output.  It also locates the taxasets that were used in the analysis.
}
\value{
A data.frame with the censored ratetest analysis results
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
\code{\link{runCensored}}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

