\name{brownie.freqs}
\alias{brownie.freqs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
List frequencies
}
\description{
List valid frequencies used in different brownie analyses
}
\usage{
brownie.freqs(with.desc = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{with.desc}{
%%     ~~Describe \code{with.desc} here~~
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (with.desc = F) 
{
    brownie_freqs_options = c("EMPIRICAL", "EQUILIBRIUM", "UNIFORM", 
        "OPTIMIZE", "SET")
    brownie_freqs_description = c("Frequencies based on empirical distribution of states at the tips", 
        "The frequencies expected with the optimized rate matrix given infinitely long branches", 
        "Frequencies all the same", "Frequencies are optimized as part of the model", 
        "User-specified frequencies")
    if (with.desc) 
        brownie_freqs_options = cbind(brownie_freqs_options, 
            brownie_freqs_description)
    return(brownie_freqs_options)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ methods }% __ONLY ONE__ keyword per line
