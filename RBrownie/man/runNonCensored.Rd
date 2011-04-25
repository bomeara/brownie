\name{runNonCensored}
\alias{runNonCensored}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add and Run NonCensored ratetest command to brownie object
}
\description{
This function runs a NonCensored ratetest on a brownie object using various models.  Commands currently in the brownie object are ignored.  
}
\usage{
runNonCensored(brobj, 
	outfile = NULL, 
	brfile = NULL, 
	usestate='simmap_state',
	models = brownie.models.continuous()[1], 
	taxsets = character(0), 
	states = character(0), 
	changes = character(0), 
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
Add a number of different continuous models. (see \code{brownie.models.continuous} for options)
}
  \item{taxsets}{
Taxasets (as character strings) to use in this NonCensored rate test.
}
  \item{states}{
This options is passed to \code{addModel}.
}
  \item{changes}{
This option is passed to \code{addModel}.
}
  \item{\dots}{
This option is passed to \code{addNonCensored}
}
}
\details{
See \code{addNonCensored} for more information about the discrete character evolution analysis itself.  This command simply combines that one with \code{run.analysis} and \code{read.continuous.output}.  

For each model one 'cont' command is add to the brownie object.  Before any 'cont' commands are added to the brownie object internally the current commands are cleared out.  This is very important behavior to take note of:  most runANALYSIS commands clear the commands (using \code{clearCommands}) before adding new ones and running the file.
}
\value{
A data.frame with the analysis output.  
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
\code{\link{addNonCensored}},
\code{\link{addModel}},
\code{\link{run.analysis}},
\code{\link{summaryCont}}
}
\examples{
# Noncensored ("cont")
\dontrun{
data(parrot)
junkrun=runNonCensored(parrot,
	brfile="cont_junk.txt",
	models=brownie.models.continuous()[1:2],
	treeloop=TRUE,
	charloop=TRUE)
summaryCont(junkrun)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

