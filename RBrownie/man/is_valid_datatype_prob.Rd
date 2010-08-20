\name{is_valid_datatype_prob}
\alias{is_valid_datatype_prob}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check for valid datatype probability
}
\description{
For future use
}
\usage{
is_valid_datatype_prob(object)
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
%%  ~~who you are~~
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
    if (object@type != discData()) 
        return(FALSE)
    if (object@ln) {
        if (object@prob > 0) 
            retval = F
    }
    else {
        if (object@prob > 1 || object < 0) 
            retval = F
    }
    return(retval)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
