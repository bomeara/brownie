\name{plot.censored}
\alias{plot.censored}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot a censored tree
}
\description{
This function takes in a brownie object and taxaset index and plots the tree, 'censoring' branches like they might be done in a censored ratetest run by brownie.  Essentially, the branch from the most recent common ancestor of the specified taxaset to it's ancestor is "censored" and, if the taxaset is paraphyletic, any internal branches within the paraphyly connecting excluded nodes.  (see example)
}
\usage{
plot.censored(x, taxind, taxcol = "red", taxlwd = 1, excol = "grey", exlwd = 1, rmlty = 3, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A brownie object
}
  \item{taxind}{
A column index for the taxaset to use.
}
  \item{taxcol}{
The color to be used for the branches in the taxaset
}
  \item{taxlwd}{
The line width to be used for the brances in the taxaset
}
  \item{excol}{
The color to be used for the branches NOT in the taxaset.  Defaults to grey
}
  \item{exlwd}{
The line width to be used for the brances NOT in the taxaset
}
  \item{rmlty}{
The line type to be used for the 'censored' branches. Defaults to dashed
}
  \item{\dots}{
Options to be passed to phyextPlot
}
}
\details{
This function is a specialization of phyextPlot.  The
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
data(parrot)
plot.censored(parrot[[1]], 2, taxlwd=2, rmlty=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

