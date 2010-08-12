\name{addSubNode}
\alias{addSubNode}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add subnode to phylo4d_ext object
}
\description{
Add a subnode to phylo4d_ext object.  This will rarely be used
}
\usage{
addSubNode(x, anc, dec, position, dataf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
is a phylo4d_ext object 
}
  \item{anc}{
Node id of the ancestor on the branch on which the subnode is to be added
}
  \item{dec}{
Node id of the descendant on the branch on which the subnode is to be added
}
  \item{position}{
The position on the branch where the subnode is to be placed.  Position is currently a 2 column matrix - this is because it will eventually be used to store confidence intervals.  For now it just the proportion of the branch starting from the ancestor.
}
  \item{dataf}{
a data.frame object which mimics the structure of the the data slot.  It contains the subnode data (which is why it should look similar to data).
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, anc, dec, position, dataf) 
{
    if (!is(x, "phylo4d_ext")) {
        warning("x is not an extended phylo4d object")
        return(x)
    }
    eind = .edge.index(x, anc, dec)
    if (is.na(eind)) 
        stop("Failure to find edge from ", anc, " to ", dec, 
            "\n")
    elen = edgeLength(x)[eind]
    if (elen == 0) 
        stop("Cannot place subnode on a zero length branch")
    if (position > elen) 
        stop("Position: ", position, ", is greater than branch length: ", 
            elen, "\n")
    newdf = getEmptyDataFrame(x)
    if (is.data.frame(dataf)) {
        if (all(names(dataf) == names(newdf))) 
            newdf = rbind(newdf, dataf)
    }
    else {
        if (ncol(newdf) == length(dataf)) {
            newdf[1, ] <- dataf
        }
        else {
            newdf[1, ] <- rep(NA, ncol(newdf))
            ndf = names(dataf)
            count = 1
            for (nm in ndf) {
                nmind = which(names(newdf) == nm)
                if (!is.null(nmind) && length(nmind) != 0) 
                  newdf[1, nmind] = dataf[count]
                count = count + 1
            }
        }
    }
    x@subnode.id = append(x@subnode.id, as.integer(nTips(x) + 
        nNodes(x) + nSubNodes(x) + 1))
    x@subnode.data = rbind(x@subnode.data, newdf)
    x@subnode.branch = rbind(x@subnode.branch, edges(x)[eind, 
        ])
    x@subnode.pos = rbind(x@subnode.pos, rep(position, 2))
    return(x)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
