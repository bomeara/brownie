\name{newlabels.v15}
\alias{newlabels.v15}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate simmap v1.5 labels
}
\description{
This is a helper function for write.simmap and is chiefly for internal use.  It generates new simmap style text labels to be written to a tree string.
}
\usage{
newlabels.v15(x, usestate, splitchar)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The phylo4d_ext object to generate simmap (1.5) labels for.
}
  \item{usestate}{
Which data column should be used written in the labels.  Note: If left as NULL (default), then all columns will be used
}
  \item{splitchar}{
Should always be a comma
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A character vector containing simmap-style taxa labels
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{newlabels.v1x}},
\code{\link{write.simmap}}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
