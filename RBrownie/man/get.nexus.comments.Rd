\name{get.nexus.comments}
\alias{get.nexus.comments}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get comments from nexus text
}
\description{
This function reads in a nexus file or uses supplied text and then returns the comments.  The comments start with [ and end with ].  This function may not work with multiline comments; don't use it expecting so.
}
\usage{
get.nexus.comments(finput, text = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{finput}{
A file which is read to text
}
  \item{text}{
Character text which can be alternatively supplied
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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

