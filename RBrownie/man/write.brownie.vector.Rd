\name{write.brownie.vector}
\alias{write.brownie.vector}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Return a vector formatted for the brownie core
}
\description{
This function can assist users in creating vectors which are compatible with the brownie core.  The inputted vector is made compatible with certain brownie commands, namely the 'statevector' option for the 'discrete' command
}
\usage{
write.brownie.vector(vect, rescale = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vect}{
The numeric vector to convert
}
  \item{rescale}{
Should the numbers in the vector be scale so that they sum to 1?
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A character string containing the formatted vector
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
vect = c(0.9,0.1) # frequencies
bvect = write.brownie.vector(vect)
data(geospiza_ext)
geotmp = addDiscrete(geospiza_ext,freq="SET",statevector=bvect )
commands(geotmp)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

