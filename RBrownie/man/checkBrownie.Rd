\name{checkBrownie}
\alias{checkBrownie}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check if brownie object is valid. 
}
\description{
Internal use only
}
\usage{
checkBrownie(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
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
function (object) 
{
    retval = TRUE
    if (length(object@commands) != 0 && !all(object@commands == 
        "")) 
        retval = retval && checkCommands(object@commands)
    if (length(object@datatypes) != 0 && !all(object@datatypes == 
        "")) 
        retval = retval && checkDataTypes(object@datatypes)
    return(retval && checkPhylo4(object))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

