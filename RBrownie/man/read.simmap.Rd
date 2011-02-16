\name{read.simmap}
\alias{read.simmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read simmap-formatted newick strings (version 1.0)
}
\description{
This function takes in a text string which contains a simmap-formatted newick string.  '{' and '}' are used in these strings to indicate subnodes and character states.
}
\usage{
read.simmap(file = "", text = NULL, vers = 1.1, as.num = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
This is not being used
}
  \item{text}{
The simmap newick string
}
  \item{vers}{
Which SIMMAP format to use.  vers=1.0 uses ';' to split up tokens within a simmap block.  vers=1.1 uses ':'.  See examples
}
  \item{as.num}{
Not currently used.
}
  \item{\dots}{
Not currently used.
}
}
\details{

}
\value{
A phylo4d_ext object
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
# SIMMAP version 1.0
rd.simmap = "((Taxon1:{A,0.1; C,0.1}, Taxon2:{T,0.1; C,0.1}):{C,0.5}, Taxon3:{C,0.4} );"
read.simmap(text=rd.simmap,vers=1.0)
# SIMMAP version 1.1
rd.simmap = "((Taxon1:{A,0.1: C,0.1}, Taxon2:{T,0.1: C,0.1}):{C,0.5}, Taxon3:{C,0.4} );"
read.simmap(text=rd.simmap,vers=1.1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

