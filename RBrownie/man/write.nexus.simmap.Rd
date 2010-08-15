\name{write.nexus.simmap}
\alias{write.nexus.simmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Write phylo4d_ext object or list to a nexus file 
}
\description{
This function is similar to \code{write.nexus} from the ape package.  If the there are subnodes in the object then these are written as SIMMAP-formatted (v1.0) strings.  
}
\usage{
write.nexus.simmap(obj, file = "", translate = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
A phylo4d_ext object or list of these objects
}
  \item{file}{
The file to output to
}
  \item{translate}{
Should taxa names be translated (not written in every newick string) or not.
}
}
\details{
This function will eventually support more subnode output formats other than just SIMMAP v1.0.  Work is in progress on this.  But for now, that is the only format that the brownie core reads and this function was written mainly to generate brownie files.

It borrows a lot of it's code and structure from the ape function \code{write.nexus}
}
\value{
Nothing
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
\code{\link{read.nexus.simmap}}
}
\examples{
\dontrun{data(geospiza_ext)}
\dontrun{write.nexus.simmap(geospiza_ext,file="geosimmap.nex")}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

