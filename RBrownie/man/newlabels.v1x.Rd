\name{newlabels.v1x}
\alias{newlabels.v1x}

\title{
Generate simmap v1.0 or v1.1 labels
}
\description{
This is a helper function for write.simmap and is chiefly for internal use.  It generates old simmap style text labels to be written to a tree string.
}
\usage{
newlabels.v1x(x, usestate, splitchar, write.nas = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The phylo4d_ext object to generate simmap (1.0 or 1.1) labels for.
}
  \item{usestate}{
Which data column should be used written in the labels.  Note: Only the first will be used
}
  \item{splitchar}{
Which character should be used to split subnodes (for 1.0 -> ';' and for 1.1 -> ':')
}
  \item{write.nas}{
Should \code{NA} values be written in the string (TRUE) or ignored (FALSE)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A character vector containing simmap-style taxa labels
}
\author{
J. Conrad Stack
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{newlabels.v15}},
\code{\link{write.simmap}}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }% __ONLY ONE__ keyword per line
