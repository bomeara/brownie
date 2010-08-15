\name{summaryRatetest}
\alias{summaryRatetest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarize the results of a Censored rate test.
}
\description{
Summarize the results of a Censored rate test.
}
\usage{
summaryRatetest(ratedf, txt = NULL, short = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ratedf}{
The data.frame returned by a Non-censored rate test
}
  \item{txt}{
Alternatively, pass in text directly from a censored rate test which this function will attempt to turn into a data.frame
}
  \item{short}{
If TRUE, then a summary of the results are returned.
}
}
\details{
The results that are printed show the support for two models: Model A uses one rate for the brownian motion parameter and Model B uses multiple rates.
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
junkrun = runCensored(parrot,taxset="intrajoint",reps=1000,charloop=TRUE)
summaryRatetest(junkrun)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

