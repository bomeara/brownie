\name{addHS}
\alias{addHS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add Heuristic search command to brownie object
}
\description{
This function adds a character string to the the commands slot of a brownie object (or list of brownie objects) which instructs the brownie program to run a heuristic search.
}
\usage{
addHS(obj, nreps = NULL, rearrlimit = NULL, maxnumspp = NULL, minnumspp = NULL, minsamp = NULL, structwt = NULL, pthreshold = NULL, subsample = NULL, movefreq = NULL, file = NULL, showtries = NULL, coal = NULL, ms = NULL, aic.mode = NULL, gridwidth = NULL, gridsize = NULL, maxrecursions = NULL, branch.export = NULL, usetempfile = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
A brownie object (or list of brownie objects)
}
  \item{nreps}{
Number of random starting species trees to use
}
  \item{rearrlimit}{
Limit search to X rearrangements for each nrep
}
  \item{maxnumspp}{
The maximum number of species to split the samples into (only relevant if SppNumFixed==No)
}
  \item{minnumspp}{
The minimum  number of samples per species (NOT the number of species)
}
  \item{minsamp}{
The minimum  number of samples per species (NOT the number of species)
}
  \item{structwt}{
StructWt
}
  \item{pthreshold}{
%%     ~~Describe \code{pthreshold} here~~
}
  \item{subsample}{
%%     ~~Describe \code{subsample} here~~
}
  \item{movefreq}{
Sets the relative proportion of times to try
	1) Species tree branch swaps
	2) Moving samples from one species to another
	3) Increasing the number of species
	4) Decreasing the number of species
	5) Attempt to reroot the species tree
If these don't sum to one, the program will automatically correct this.
}
  \item{file}{
The brownie core write the output of the analysis to this file (optional)
}
  \item{showtries}{
%%     ~~Describe \code{showtries} here~~
}
  \item{coal}{
Not currently used
}
  \item{ms}{
Not currently used
}
  \item{aic.mode}{
Not currently used
}
  \item{gridwidth}{
affects grid search
}
  \item{gridsize}{
affects grid search
}
  \item{maxrecursions}{
affects grid search
}
  \item{branch.export}{
if set to 0, returns a single estimate of the best branch lengths on the species tree. If set to 1, returns a table of equally-good branch lengths. If set to 2, returns a table of the best and neighboring branch lengths, suitable for doing a contour plot of score versus branch lengths
}
  \item{usetempfile}{
Should a tempfile be written to (instead of file).  Not currently used.
}
}
\details{
Do a heuristic search.

From brownie command-line:
Keyword ---- Option type --------------------------- Current setting --
NReps        <integer-value>                         5
RearrLimit   <integer-value>|None                    None
MaxNumSpp    <integer-value>                         100
MinNumSpp    <integer-value>                         1
MinSamp      <integer-value>                         3
StructWt     <double>                                0.5
PThreshold   <double>                                1.0
Subsample    <double>                                2.0
MoveFreq     (number number number number number [number])    (0.8 0.01 0.01 0.01 0.17 0.0 )
File         <file-name>                             besttrees.tre
ShowTries    No|Yes                                  No
COAL:        No|Yes                                  No
MS:          No|Yes                                  No
AIC_mode      0|1|2|3|4                              4
GridWidth     <double>                               10.0
GridSize      <integer-value>                        15.0
MaxRecursions <integer-value>                        20
Branch_export 0|1|2                                  2
}
\value{
Returns the original object (obj) with ratetest command added.
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
ptmp = clearCommands(parrot)
ptmp = addHS(ptmp,file="besttrees.tre2")
commands(ptmp)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
