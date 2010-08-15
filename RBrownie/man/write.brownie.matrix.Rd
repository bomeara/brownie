\name{write.brownie.matrix}
\alias{write.brownie.matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Return a matrix formatted for the brownie core
}
\description{
This function can assist users in creating matrices which are compatible with the brownie core.  The inputted matrix is essentially flattened, make it compatible with certain brownie commands, namely 'discrete'
}
\usage{
write.brownie.matrix(mat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{
%%     ~~Describe \code{mat} here~~
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
\code{\link{addDiscrete}}
}
\examples{
testmat = matrix(c("","a",
			"b",""), byrow=T,nrow=2)
bmat = write.brownie.matrix(testmat)
data(geospiza_ext)
geotmp = addDiscrete(geospiza_ext,model="USER",ratemat=bmat)
commands(geotmp)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

