\name{addNonCensored}
\alias{addNonCensored}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add Noncensored rate test command to brownie object
}
\description{
This function adds a character string to the the \code{commands} slot of a brownie object (or list of brownie objects).  This character string instructs the brownie program to run a non-censored ratetest.
}
\usage{
addNonCensored(obj, 
	file = NULL, 
	taxset = character(0), 
	treeloop = FALSE, 
	charloop = FALSE, 
	fileappend = FALSE, 
	filereplace = FALSE, 
	usetempfile = FALSE, 
	...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
A brownie object (or list of brownie objects)
}
  \item{file}{
The brownie core write the output of the analysis to this file (optional)}
  \item{taxset}{
Specifies which taxa (see \code{taxasets(obj)}) in the tree should be used in the analysis, usually when running a non-censored test it is all the species in the tree but it can be used to remove species with missing data.
}
  \item{treeloop}{
Should the analysis loop across the trees?  To use a specific tree then run the \code{addChoose(tree=<TREENUM>)} before this command.
}
  \item{charloop}{
Should the analysis loop across characters?  To use a specific character then run the \code{addChoose(char=<CHARNUM>)} before this command.
}
  \item{fileappend}{
If the output file exists, appends to it rather than overwrites it. Turned on by default.
}
  \item{filereplace}{
If set to yes, if the output file already exists it will be quietly replaced.
}
  \item{usetempfile}{
Should a tempfile be written to (instead of file).  Not currently used.
}
  \item{\dots}{
Not currently used.
}
}
\details{
This is still a bit untested, but has been used successfully (and has been tested somewhat). Basically, parameters can be "painted" onto the branches of the tree. For example, you could paint most of the tree with one rate parameter, and just one clade (including half the branch leading up to it) with a different rate parameter. This \dQuote{painting} can be done in one of two ways: 1) loading a file in SIMMAP (Bollback, 2006) format with the numbers 0-9 mapped on the tree, or 2) reconstructing states of a discrete trait down the tree (next section). Unlike programs like OUCH (Butler and King, 2004), which paints Ornstein-Uhlenbeck regimes on entire branches, Brownie can handle branches where the first 24% has one regime, the next 15% has another regime, and the remaining part has a third regime. The painted parameters can be different Brownian motion rate parameters or different OU means (more models are coming soon).
}
\value{
Returns the original object (obj) with CONTinuous command added.
}
\references{
From tutorial written by Samantha Price and Peter Wainwright. Presented at Bodega Bay Phylogenetics Workshop March 6-13th 2010. 
}
\author{
J. Conrad Stack
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{plot.taxaset}}
}
\examples{
data(parrot)
ptmp = clearCommands(parrot)
ptmp <- addModel(ptmp,"BM1")
ptmp <- addNonCensored(ptmp,file="parrotfish_feeingmodetest.txt",treeloop=TRUE,charloop=TRUE,taxset="all",filereplace=TRUE)
ptmp <- addModel(ptmp,"BMS")
ptmp <- addNonCensored(ptmp,file="parrotfish_feeingmodetest.txt",treeloop=TRUE,charloop=TRUE,taxset="all",fileappend=TRUE,filereplace=FALSE)
commands(ptmp)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
