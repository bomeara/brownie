\name{read.nexus.block}
\alias{read.nexus.block}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find a block in a nexus file
}
\description{
Scan a nexus file for a certain block.  This is mainly an internal function, but might be used to diagnose problems with nexus formatted files.  It doesn't attempt to process the nexus file, just searches the rawtext for blocks
}
\usage{
read.nexus.block(finput, txt = NULL, block, rm.comments = F, silent = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{finput}{
A nexus file
}
  \item{txt}{
Alternatively, look through this text
}
  \item{block}{
The block name to look for
}
  \item{rm.comments}{
Logical. Should comments be removed?  This option is experimental at the moment and does not work for multiline comments.
}
  \item{silent}{
Verbose output?
}
}
\details{
Eventually, multiline comments will be supported once the read and scanning of nexus files becomes more sophisticated.
}
\value{
A character vector of the text in the \code{block} block.
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (finput, txt = NULL, block, rm.comments = F, silent = F) 
{
    if (!is.null(txt)) {
        rawtext = txt
    }
    else {
        if (!file.exists(finput)) 
            stop("Assuming finput is a file and could not find it")
        rawtext = scan(finput, what = character(0), strip.white = T, 
            sep = "\n", quiet = silent)
    }
    inds = .get.nexus.block.inds(blockname = block, text = rawtext)
    if (length(inds) == 0) {
        if (!silent) 
            warning(paste("This file has", block, "no block"))
        return(character(0))
    }
    rawtext = rawtext[(inds[1] + 1):(inds[2] - 1)]
    if (rm.comments) {
        comment.starts = grep("^\\[", rawtext)
        comment.ending = grep("\\]$", rawtext)
        if (length(comment.starts) == length(comment.ending)) {
            comment.pairs = cbind(comment.starts, comment.ending)
        }
        else {
            comment.pairs = cbind(character(0), character(0))
        }
        if (nrow(comment.pairs) > 0) {
            for (pair in seq(nrow(comment.pairs))) {
                print(pair)
                rawtext = rawtext[-(comment.pairs[pair, 1]:comment.pairs[pair, 
                  2])]
            }
        }
    }
    return(rawtext)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

