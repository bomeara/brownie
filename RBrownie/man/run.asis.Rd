\name{run.asis}
\alias{run.asis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Execute a brownie object
}
\description{
This function is the sister function to run.analysis.  Instead of running a file directly it writes \code{brobj} to either \code{brfile} or a temporary file and then executes.  
}
\usage{
run.asis(brobj, brfile = NULL, logfile = NULL, debugging = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{brobj}{
A brownie object or list of brownie objects
}
  \item{brfile}{
A file to write the brownie object out to.  If this is not set then \code{brobj} is written to a temporary file.
}
  \item{logfile}{
Not currently used
}
  \item{debugging}{
Not currently used
}
}
\details{
This function returns very raw text from the brownie core and would be good for debugging or if you wanted more control over how to execute your brownie objects.  There are 'pre-packaged' runTEST commands which make it easier to run RBRownie analyses, but used in combination with \code{addLiteral}, this command give more flexibility in how to setup and execute brownie analyses.
}
\value{
\item{treesout}{Newick text strings which the brownie core generated.  'discrete' when used with the option reconstruct=yes will generate these trees.  Often they are SIMMAP-formatted (v1.0)}
\item{textout}{Raw text of the analysis output.  Try using scan.textout and read.analysis.output to get this information into a more tractable data structure}
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
\code{\link{run.analysis}},
\code{\link{scan.textout}},
\code{\link{scan.textout}}
}
\examples{
data(parrot)
junkrun=run.asis(parrot)
cuttxt = scan.textout(junkrun$textout[1]) # post-process the text
read.analysis.output(txt=cuttxt) # read output into a data.frame
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

