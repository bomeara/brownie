\name{write.simmap.new}
\alias{write.simmap.new}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Wrapper for write.simmap
}
\description{
Just a wrapper function for write.simmap which instructs write.simmap to use simmap version 1.5
}
\usage{
write.simmap.new(x, usestate = NULL, file = "", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A tree object compatible with write.simmap
}
  \item{usestate}{
Which state to output?  If left as NULL, then all the data columns are used
}
  \item{file}{
Where out output the simmap style tree string to?
}
  \item{\dots}{
Other arguments which are passed on to write.simmap
}
}

\author{
J. Conrad Stack
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{write.simmap}}
}

\keyword{ methods }

