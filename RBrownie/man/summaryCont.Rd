\name{summaryCont}
\alias{summaryCont}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarize the results of a Non-censored rate test.
}
\description{
This prints out the data.frame returned from a Non-censored rate test and is very helpful when attempting to understand the results.
}
\usage{
summaryCont(contdf, txt = NULL, short = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{contdf}{
The data.frame returned by a Non-censored rate test
}
  \item{txt}{
Alternatively, pass in text directly from a Non-censored rate test which this function will attempt to turn into a data.frame
}
  \item{short}{
If TRUE, then a summary of the results are returned.
}
}
\details{
The results that are printed show to different model comparisons (AIC and AICc values)
}
\value{
Print out to screen
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
junkrun=runNonCensored(parrot,brfile="cont_junk.txt",models=brownie.models.continuous()[1:2],treeloop=TRUE,charloop=TRUE)
summaryCont(junkrun)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

