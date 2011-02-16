\name{plot.taxaset}
\alias{plot.taxaset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot a taxaset 
}
\description{
Plot a taxaset onto the tree.  Basically, this function calls phyextPlot using a column of taxaset binary data as the usestate and excol and taxcol as the state.cols. 
}
\usage{
plot.taxaset(x, taxind, taxcol = "red", taxlwd = 1, excol = "grey", exlwd = 1, blankit = F, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A brownie object
}
  \item{taxind}{
The taxaset column index to plot
}
  \item{taxcol}{
The color to use for the branches in the taxaset
}
  \item{taxlwd}{
The line width to use for the branches in the taxaset
}
  \item{excol}{
The color to use for the branches NOT in the taxaset
}
  \item{exlwd}{
The line width to use for the branches NOT in the taxaset
}
  \item{blankit}{
Logical.  Should the branch connecting the taxaset to the rest of the tree be removed (not shown)
}
  \item{\dots}{
These arguments are passed to phyextPlot
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
data(parrot)
plot.taxaset(parrot[[1]], 2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

