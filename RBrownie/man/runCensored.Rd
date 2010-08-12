\name{runCensored}
\alias{runCensored}
\alias{Censored Ratetest}
\title{
Add and Run Censored ratetest command to brownie object
}
\description{
This function runs a Censored ratetest on a brownie object.  Commands currently in the brownie object are ignored.  
}
\usage{
runCensored(brobj, brfile = NULL, file = NULL, ...)
}
\arguments{
  \item{brobj}{
A brownie object (or list of brownie objects)
}
  \item{brfile}{
A file to write brobj to (before running) (optional).
}
  \item{file}{
The brownie core write the output of the ratetest to this file (optional)
}
  \item{\dots}{
Arguments which are passed to addCensored
}
}
\details{
The censored ratetest compares the rates of evolution in a monophyletic clade to another monophyletic of paraphyletic clade (i.e. rate change only occurs at nodes).  It is useful for asking questions about whether a key character has lead to increased rates of evolution.  It is 'censored' as it removes the branches that join the two clades together as it is not clear which clade to assign these branches to.
}
\value{
A data.frame with the results of the model run.  The data.frame may have different columns depending on the options chosen, but it will contain a row for each tree and character used in the analysis.  It also includes the maximum likelihood values of the brownian motion rates inferred for each set of taxa.  The MLE of the state of the root node of each clade and the ML value itself (as log-likelihood) is also return.  It also returns the relative support (using different criteria) for one of two rate models.  Model A is a 1-rate model. Model B is a 2-rate model.

	The 'sigcodes' variable gives a summary of the support for models A and B.  RBrownie uses the notation A,B,a,b to indicate strength of support for a single rate ('A's) or two rate ('B's) model, weak ('a') or strong ('A') support is indicated by capitalization.
}

\references{
%% ~put references to the literature/web site here ~
}
\author{
J.Conrad Stack
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{summaryRatetest}},
\code{\link{plot.ratetest}}
}
\examples{
# run Censored ratetest
data(parrot)
brownierun = runCensored(parrot,taxset="intrajoint",reps=1000,charloop=T)
summaryRatetest(brownierun)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
