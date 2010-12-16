\name{areTaxaMutex}
\alias{areTaxaMutex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Are taxasets mutually exclusive
}
\description{
This function determines whether or not an indefinitate number of taxasets are mutually exclusive.  To put it another way, do any of the taxasets specified overlap?
}
\usage{
areTaxaMutex(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A brownie object
}
  \item{\dots}{
An indefinite number of taxaset column indices.  See example below.
}
}
\details{
This function is still relatively experimental.
}
\value{
TRUE if the taxaset column indices in \code{...} are mutually exclusive.
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

data(parrot)
brobj = parrot[[10]]
newtaxaset = c("Sparisoma_radians",
"Sparisoma_chrysopterum",
"Sparisoma_viride",
"Sparisoma_aurofrenatum",
"Sparisoma_rubripinne",
"Sparisoma_atomarium")

# This 'replace' function actually ADDS a taxaset
taxasets(brobj,taxnames="newset") <- newtaxaset
taxasets(brobj)
areTaxaMutex(brobj,1,3) # FALSE
areTaxaMutex(brobj,2,3) # TRUE

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
