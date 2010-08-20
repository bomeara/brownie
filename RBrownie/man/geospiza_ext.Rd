\name{geospiza_ext}
\alias{geospiza_ext}
\docType{data}
\title{
Data from Darwin's finches
}
\description{
Phylogenetic tree and morphological data for Darwin's finches, in different formats
}
\usage{data(geospiza_ext)}
\format{
  The format is:
Formal class 'brownie' [package "RBrownie"] with 15 slots
}
\details{
Phylo4d object stolen from phylobase and extended into a brownie object.
}
\source{
phylobase
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(geospiza_ext)
plot(geospiza_ext)
commands(geospiza_ext)
geomod <- removeCommands(geospiza_ext,"log")
commands(geomod)
# run the commands currently in the commands block
\dontrun{run.asis(geomod)}
}
\keyword{datasets}
