\name{get.tree.weights}
\alias{get.tree.weights}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get tree weights
}
\description{
Read the tree weights from a nexus file or text.  It looks within the comments before each tree string for the weights.
}
\usage{
get.tree.weights(finput, 
	text = NULL, 
	starters = c("&W", "&lnP"), 
	splitchar = "( |=)")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{finput}{
A file to be read in
}
  \item{text}{
Alternatively, send it text itself.
}
  \item{starters}{
The character pattern which indicates the start of the weight (within a comment)
}
  \item{splitchar}{
The character pattern which separates the 'starters' from the weights themselves
}
}
\details{

}
\value{
Returns a numeric vector of the weights (if any are avialable).
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

