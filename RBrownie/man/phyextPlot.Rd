\name{phyextPlot}
\alias{phyextPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot phylo4d_ext object
}
\description{
This function plots a phylo4d_ext object.  Branches and nodes in different states can be colored.  If subnodes exist it also plots these on the tree, showing them on the exact position they actually occur.   This plotting function is an extension to treePlot in phylobase.  It was created to plot discrete data states on top of a tree.  
}
\usage{
phyextPlot(x, states, states.col, states.na = "none", datapart = 1, plot.subnodes = T, plot.points = T, line.widths, line.types, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The phylo4d_ext object
}
  \item{states}{
A vector containing all the discrete states to show on the plot.  
}
  \item{states.col}{
The colors for each state. The branches are painted these colors
}
  \item{states.na}{
If \code{states} is missing, then this is used to represent the NAs in the dataset.  It is mainly for internal use.
}
  \item{datapart}{
This is the column index for a discrete data column whose states will be painted onto the tree.  The data column specified should contain data for every node in the tree.  For any missing data the \code{states.na} will be used.
}
  \item{plot.subnodes}{
Logical specifying whether or not subnodes should be plotted (if they exist).
}
  \item{plot.points}{
Logical specifying whether or not points should be plotted at each internal node and subnode.  This can help indicate where subnodes exist on a branch.
}
  \item{line.widths}{
Integer. Line widths for each state.  Branches in this state will use this line width.
}
  \item{line.types}{
Integer. Line types for each state.  Branches in this state will use this line type.
}
  \item{\dots}{
This is passed to the phylobase function treePlot, which is the function that this function uses to place the tree on the canvas.
}
}
\details{
This is nearly an internal function.  It basically uses treePlot (to which you can pass options via ...) to plot a phylogenetic tree (with no continuous data).  It then plots discrete data on top of this plot based on the which state each branch is in.  The branch is colored based on the state of the descendant node.  

Note that treePlot is called first and actually drawn on the canvas before the colored branches are plotted on top.  This could lead to some odd plots when the using different line.widths and line.types.  To neutralize this effect try passing colors to treePlot through ... which paints the original plot the same color as the background color.

This is the default function for phylo4d_ext.
}
\value{
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

# default plot
phyextPlot(parrot[[1]])
phyextPlot(parrot[[1]],datapart=1) # same plot as the last one

# the states are 0,1 in tdata(parrot[[1]])[,1]
# use black and red as the colors.
plot(parrot[[1]],states=c(0,1),states.col=c(1,2))
# show without points
plot(parrot[[1]],states=c(0,1),states.col=c(1,2),plot.points=FALSE)
#
# Pass the option \code{show.tip.label} to treePlot
plot(parrot[[1]],states=c(0,1),states.col=c(1,2),plot.points=FALSE,show.tip.label=FALSE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

