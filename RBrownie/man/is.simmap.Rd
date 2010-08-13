\name{is.simmap}
\alias{is.simmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check SIMMAP v1.0
}
\description{
Check if text of file contains a SIMMAP-formatted (v1.0) newick file
}
\usage{
is.simmap(finput = "", text = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{finput}{
A newick file
}
  \item{text}{
Newick file text
}
}
\details{
This function looks for '{' and '}' characters in the text which would indicate that it is a SIMMAP-formatted file
}
\value{
Logical
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
phytext = "((Taxon1:{A,0.1; C,0.1}, Taxon2:{T,0.1; C,0.1}):{C,0.5}, Taxon3:{C,0.4} );"
is.simmap(text=phytext) # TRUE
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

