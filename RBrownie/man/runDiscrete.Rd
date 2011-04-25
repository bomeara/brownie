\name{runDiscrete}
\alias{runDiscrete}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add and Run a discrete command to brownie object
}
\description{
This function runs a discrete character evolution analysis on a brownie object using various models and frequencies.  Commands currently in the brownie object are ignored.  
}
\usage{
runDiscrete(brobj, 
	outfile = NULL, 
	brfile = NULL,
	usestate='simmap_state', 
	models = brownie.models.discrete()[1], 
	freqs = brownie.freqs()[1], 
	...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{brobj}{
A brownie object (or list of brownie objects)
}
  \item{outfile}{
A file to write the results of the test to.  If not specified a temporary file is used
}
  \item{brfile}{
A file to write brobj to (before running) (optional).
}
  \item{usestate}{
	This is passed to writeBrownie
}
  \item{models}{
Add a number of different discrete models. (see \code{brownie.models.discrete} for options)
}
  \item{freqs}{
Add a number of different discrete state frequency models.  (see \code{brownie.freqs} for more options)
}
  \item{\dots}{
All these arguments are passed to \code{addDiscrete}.
}
}
\details{
See \code{addDiscrete} for more information about the discrete character evolution analysis itself.  This command simply combines that one with \code{run.analysis} and \code{read.discrete.output}.  

For each model / freq combination one discrete command is add to the brownie object.  Before any discrete commands are added to the brownie object internally the current commands are cleared out.  This is very important behavior to take note of:  most runANALYSIS commands clear the commands (using \code{clearCommands}) before adding new ones and running the file.
}
\value{
%%  If it is a LIST, use
\item{stats}{The statistical output of the run}
\item{trees}{The reconstructed trees (if reconstruct=T) for each model/freq combination}
\item{cmds}{The commands, in brownie command-line format, which were used}

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
\code{\link{addDiscrete}},
\code{\link{brownie.freqs}},
\code{\link{run.analysis}},
}
\examples{
\dontrun{
data(geospiza_ext)
junkrun=runDiscrete(geospiza_ext,
	models=c("nonrev","nonrev"),
	freqs=c("unif","equilib"),
	reconstruct=TRUE)
plot(junkrun$trees[[1]]) # plot one of the reconstructed trees
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

