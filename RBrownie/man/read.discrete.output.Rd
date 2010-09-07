\name{read.discrete.output}
\alias{read.discrete.output}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read brownie output from discrete analysis
}
\description{
This function processes output returned from the brownie core before passing it on to read.analysis.output.
}
\usage{
read.discrete.output(filename, txt = NULL, show.warn = FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
A file with brownie output
}
  \item{txt}{
Text containing brownie output
}
  \item{show.warn}{
Should internal warnings be shown?  If you are using trees with a small number of taxa or characters then this should be turned on.
}
  \item{\dots}{
Not currently used
}
}
\details{
This function removes WARNINGS returned by the brownie core output before passing the rest of it on to read.analysis.output.  This is also somewhat of an internal function
}
\value{
\item{stats}{A data.frame of the results of the analysis}
\item{trees}{The trees that might be returned by the analysis}
\item{cmds}{The discrete commands as character strings that were used. }
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
\code{\link{read.analysis.output}},
\code{\link{runDiscrete}}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

