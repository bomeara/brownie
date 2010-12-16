\name{read.analysis.output}
\alias{read.analysis.output}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read brownie output to data.frame
}
\description{
Method to scan a brownie program's output and return it as a data frame. Brownie generally output multiple header lines along with the actual raw numbers, these need to be removed before the data frame is returned.  This is chiefly an internal function.
}
\usage{
read.analysis.output(filename, txt = NULL, rowsep = "\n", colsep = "\t")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
The filename which contains the brownie analysis output 
}
  \item{txt}{
The character vector which contains the brownie analysis output.  It is an alternative to filename
}
  \item{rowsep}{
This character or pattern is what the function uses to distinguish tokens within a row
}
  \item{colsep}{
This character or pattern is what the function uses to distinguish rows
}
}
\details{
This functions parses text rows of brownie analyses looking for data lines and header lines.  There would be multiple header lines in brownie output if say multiple models were tested. 

Some analyses return different headers.  The final data frame contains the union of all headers.  If a certain model does not contain information for a header then NAs are put in as place holders. Again, this is chiefly and internal function.
}
\value{
A data.frame object
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

\seealso{

}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

