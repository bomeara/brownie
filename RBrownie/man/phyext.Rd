\name{phyext-methods}
\docType{methods}
\alias{phyext}
\alias{phyext-methods}
\alias{phyext,character-method}
\alias{phyext,list-method}
\alias{phyext,phylo-method}
\alias{phyext,phylo4-method}
\alias{phyext,phylo4d-method}
\alias{phyext,phylo4d_ext-method}
\title{A function which returns a phylo4d_ext object}
\description{
 This is a shorthand way of creating phylo4d_ext objects
}
\section{Methods}{
\describe{

\item{\code{signature(x = "character")}}{ This is not implemented yet. }

\item{\code{signature(x = "list")}}{ Convert a list of phylo,phylo4,or phylo4d objects into a list of phylo4d_ext objects }

\item{\code{signature(x = "phylo")}}{ Convert to phylo4d_ext }

\item{\code{signature(x = "phylo4")}}{ Convert to phylo4d_ext }

\item{\code{signature(x = "phylo4d")}}{ Convert to phylo4d_ext }

\item{\code{signature(x = "phylo4d_ext")}}{ Returns the argument without modification }
}}
\keyword{methods}

