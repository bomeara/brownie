\name{read.simmap.new}
\alias{read.simmap.new}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read simmap-formatted newick strings (version 1.5)
}
\description{
This function takes in a text string which contains a simmap-formatted newick string.  '[' and ']' are used in these strings to indicate subnodes and character states.
}
\usage{
read.simmap.new(file = "", text = NULL, specialpatt = character(0))
}

\arguments{
  \item{file}{
A file containing a newick string with SIMMAP 1.5 style comments
}
  \item{text}{
The simmap newick string
}
  \item{specialpatt}{
A pattern to match against variable names.  NOTE: For now, it should almost always be set to '.*' and needs to be further tested before anything else should be tried.
}

}
\details{
This is still a work in progress.  
}
\value{
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
treestr = "(('Taxon1':[&simmap_state={A,0.2,A,0.5,C}]0.8,'Taxon2':[&simmap_state={T,0.2,C}]0.35)Internal0000001:[&simmap_state={C}]0.5,Taxon3:[&simmap_state={C}]0.4)Internal0000002:[&simmap_state={G}]0;"
tr = read.simmap.new(text=treestr, specialpatt=".*")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

