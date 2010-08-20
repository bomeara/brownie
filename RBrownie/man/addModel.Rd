\name{addModel}
\alias{addModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add model selection command to brownie object
}
\description{
This function adds a character string to the the commands slot of a brownie object (or list of brownie objects) which instructs the brownie program to set a certain model for further analysis.
}
\usage{
addModel(obj, type = brownie.models.continuous()[1], states = NULL, changes = NULL, model.state = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
A brownie object (or list of brownie objects)
}
  \item{type}{
A model of the following type:
BM1 = Brownian motion, one rate parameter
BMS = Brownian motion, with different rate parameters for each state on a tree
BMC = Brownian motion, with one rate parameter for branches with state changes and another for branches without changes
OUSM = Ornstein-Uhlenbeck with one mean per discrete state (attraction and rate parameters constant across tree)
OUCM = Ornstein-Uhlenbeck with independent means on branches with and without changes in a discrete character (attraction and rate parameters constant across tree)
(see brownie.models.continuous() for more information)
}
  \item{states}{
State vector allows restrictions, so that character states 0 and 2, for example, may be viewed by the program as identical. Example: state=(0 1 0 2 3 4 5 6 7 8)
}
  \item{changes}{
Changes is the maximum number of times a particular character state can be present on a root to tip lineage. For example, if a taxon sister to all other taxa starts in state 0, changes to state 1, and then changes to state 0, state 0 has beeen present on that branch twice.
}
  \item{model.state}{
Which data.frame value should you use.  The default is to use the column name headed 'simmap_state'
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
function (obj, type = brownie.models.continuous()[1], states = NULL, 
    changes = NULL, model.state = NULL) 
{
    if (!is.list(obj)) {
        obj = list(obj)
    }
    newcmd = list(command = NULL, options = matrix(NA, ncol = 2, 
        nrow = 0))
    newcmd$command = "model"
    newcmd = addcmd.model.continuous(newcmd, type)
    if (!is.null(states)) {
        dat = NULL
        if (is.null(model.state)) {
            colind = head(which(datatypes(obj[[1]]) == discData()), 
                1)
            if (length(colind) == 0) 
                stop("Could not find data column to use.")
            dat = tdata(obj[[1]], "tip")[, colind, drop = T]
        }
        else {
            if (is.integer(model.state)) {
                colind = model.state
            }
            else {
                colind = which(colnames(tdata(obj[[1]], "tip")) == 
                  model.state)
            }
            dat = tdata(obj[[1]], "tip")[, colind, drop = T]
        }
        newcmd = addcmd.states(newcmd, states, dat)
    }
    if (!is.null(changes) && !is.na(as.integer(changes)) && as.integer(changes) > 
        0) 
        newcmd = addcmd.literal(newcmd, "changes", as.character(changes))
    for (ii in seq(length(obj))) {
        commands(obj[[ii]], add = T) <- write.brownie.string(newcmd)
    }
    if (length(obj) == 1) 
        obj = obj[[1]]
    return(obj)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

