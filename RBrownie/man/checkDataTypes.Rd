\name{checkDataTypes}
\alias{checkDataTypes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check datatypes
}
\description{
Check if brownie datatypes are valid
}
\usage{
checkDataTypes(dtypes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dtypes}{
%%     ~~Describe \code{dtypes} here~~
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
function (dtypes) 
{
    return(all(dtypes \%in\% brownie.datatypes()))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

