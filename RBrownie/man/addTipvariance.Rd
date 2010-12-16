\name{addTipvariance}
\alias{addTipvariance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add tipvariance command
}
\description{
This function adds a character string to the the commands slot of a brownie object (or list of brownie objects) which instructs the brownie program use specific tip variances.
}
\usage{
addTipvariance(obj, type = brownie.tvtypes()[1])
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
A brownie object (or list of brownie objects)
}
  \item{type}{
Observations of taxon means have uncertainty. By default, the program assumes that these values are known exactly. If you have measured these variances for each taxon (so that, for each measured character, the data matrix has the mean followed by the variance), choose Given. If you want estimate one tip variance across all the taxa, choose Same. If you want to assume no tip variance (one fewer parameter), choose None.

(see brownie.tvtypes() for more information)
}
}
\details{
Note that if the tip variance type is set to given, values in the matrix are assumed to consist of a taxon average followed by the corresponding variance, so the matrix might be:
\tabular{llllllll}{
taxon1 \tab  char1_mean \tab char1_variance \tab char2_mean \tab char2_variance \tab char3_mean \tab char3_variance \tab ... \cr
taxon2 \tab char1_mean \tab char1_variance \tab char2_mean \tab char2_variance \tab char3_mean \tab char3_variance \tab ... \cr
taxon3 \tab char1_mean \tab char1_variance \tab char2_mean \tab char2_variance \tab char3_mean \tab char3_variance \tab ... \cr
}
This differs somewhat from the NEXUS specification (Maddison et al, 1997), which would use an ITEMS=(AVERAGE VARIANCE) command and a somewhat different matrix.
}
\value{
Returns the original object (\code{obj}) with 'tv' command added to its \code{commands} slot.
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
ptmp = addTipvariance(ptmp,"given")
commands(ptmp)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
