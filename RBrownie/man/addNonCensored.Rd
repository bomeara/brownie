\name{addNonCensored}
\alias{addNonCensored}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add Noncensored rate test command to brownie object
}
\description{
This function adds a character string to the the commands slot of a brownie object (or list of brownie objects) which instructs the brownie program to run a non-censored ratetest.
}
\usage{
addNonCensored(obj, file = NULL, taxset = character(0), treeloop = FALSE, charloop = FALSE, fileappend = FALSE, filereplace = FALSE, usetempfile = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
A brownie object (or list of brownie objects)
}
  \item{file}{
The brownie core write the output of the analysis to this file (optional)}
  \item{taxset}{
Specifies which taxa (from taxasets(obj)) in the tree should be used in the analysis, usually when running a non-censored test it is all the species in the tree but it can be used to remove species with missing data.
}
  \item{treeloop}{
Should the analysis loop across the trees?  To use a specific tree then run the addChoose(tree=<TREENUM>) before this command.
}
  \item{charloop}{
Should the analysis loop across characters?  To use a specific character then run the addChoose(char=<CHARNUM>) before this command.
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
Brownie uses simmap (v1.0) output format to read in ancestral character state reconstrutions.  These reconstructions can be used with the Non-censored rate test command as a means to address questions about whether certain character states influence the evolutionary rates (under different models of evolution)
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
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
