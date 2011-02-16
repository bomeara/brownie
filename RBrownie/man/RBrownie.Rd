\name{RBrownie}
\alias{RBrownie}
\alias{RBrownie}
\docType{package}
\title{
Continuous and discrete ancestral character reconstruction and evolutionary rate tests.
}
\description{
This R package wraps the Brownie c++ library, which is a phylogenetic library which performs evolutionary rate tests and ancestral character reconstructions.
}
\details{
\tabular{ll}{
Package: \tab RBrownie\cr
Type: \tab Package\cr
Version: \tab 0.1.0\cr
Date: \tab 2010-09-10\cr
URL: \tab http://code.google.com/p/brownie/, http://www.conradstack.info\cr
License: \tab GPL-2\cr
Collate: \tab processNexus.R class-phylo4d_ext.R methods-phylo4d_ext.R setAs-phyext.R phyextPlots.R ci.R processBrownie.R class-brownie.R methods-brownie.R setAs-brownie.R browniePlots.R readWrite.R rcpp_hello_world.R Rpack-internal.R zzz.R constructAnalyses.R readAnalyses.R class-datatypes.R\cr
LazyLoad: \tab yes\cr
Depends: \tab methods, Rcpp (>= 0.8.2), phylobase (>= 0.6.2), ape (>= 2.5)\cr
LinkingTo: \tab Rcpp\cr
SystemRequirements: \tab GNU make, GNU GSL (>= 1.10)\cr
Packaged: \tab 2011-02-16 19:35:56 UTC; root\cr
Built: \tab R 2.12.0; i486-pc-linux-gnu; 2011-02-16 19:36:35 UTC; unix\cr
}

Index:
\preformatted{
addCensored             Add Censored ratetest to brownie object
addChoose               Add choose command to a brownie object (or list
                        of brownie objects)
addDiscrete             Add Discrete character models and
                        reconstructions to brownie object
addEndLog               Add log stop command to brownie object
addHS                   Add Heuristic search command to brownie object
addLiteral              Add command string to brownie object
addModel                Add model selection command to brownie object
addNonCensored          Add Noncensored rate test command to brownie
                        object
addNumopt               Add numerical optimization command to brownie
                        object
addOpt                  Add Noncensored rate test command to brownie
                        object
addStartLog             Add log start command to brownie object
addSubNode              Add subnode to phylo4d_ext object
addTipvariance          Add tipvariance command
addcmd.binary           Adding to commands block
addcmd.freq             Adding to commands block
addcmd.literal          Adding to commands block
addcmd.model.continuous
                        Adding to commands block
addcmd.model.discrete   Adding to commands block
addcmd.states           Adding to commands block
addcmd.taxaset          Adding to commands block
addcmd.tvtype           Adding to commands block
areTaxaMono             Do taxa form a monophyletic clade
areTaxaMutex            Are taxasets mutually exclusive
areTaxaPara             Do taxa form a paraphyletic clade
as.contData             Convert to continuous data
as.discData             Convert to discrete data
brownie-class           Class "brownie"
brownie-methods         Methods for generic brownie constructor
                        function
brownie.datatypes       List datatypes
brownie.freqs           List frequencies
brownie.models.continuous
                        List valid continuous models
brownie.models.discrete
                        List valid discrete models
brownie.tvtypes         List tipvariance types
check.empty             Check character vectors
check.freq              Check frequency
check.statevector       Check statevector
checkBrownie            Check if brownie object is valid.
checkCommands           Check commands in a brownie object
checkDataTypes          Check datatypes
checkMono               Check for monophyly
checkPara               Check for paraphyly
checkval.dummy          Check value
checkval.integer        Check value
checkval.model.continuous
                        Check value
checkval.model.discrete
                        Check value
checkval.ratemat        Check value
checkval.yesno          Check value
clearCommands-methods   Clear commands
collapse.singletons     Collapse singleton nodes
collapse.to.singles     Collapse zero-length branches
contData                Return datatype string representing the
                        continuous datatype
discData                Return datatype string representing the
                        discrete datatype
execute.brownie         Execute a brownie-formatted nexus file
expand.singles          Expand singletons
genericData             Return datatype string representing the generic
                        datatype
geospiza_ext            Data from Darwin's finches
get.left.right          Get right and left tokens
get.nexus.comments      Get comments from nexus text
get.nodenames           Read node names from newick text
get.tree.weights        Get tree weights
getEmptyDataFrame       Get empty data.frame
getSubNodeData          Return subnode data
getSubNodeEdgeInds      Return subnode edge indices
getSubNodePosish        Return the subnode position
has.block               Check for nexus block
has.characters2         Check for CHARACTERS2 block
has.weights             Check for tree weights
hasCommands             Check for commands
hasTaxasets-methods     Check if there are taxasets
haschar                 Does character exist
is.binary               Is data binary
is.contData             Is data continuous
is.discData             Is data discrete
is.simmap               Check SIMMAP v1.0
is.taxaname.internal    Check for internal taxanames
is_valid_datatype       Check for valid datatype
is_valid_datatype_bound
                        Check for valid datatype limits
is_valid_datatype_prob
                        Check for valid datatype probability
isback                  isback
match.order             Matching orders
nSubNodes               Get number of subnodes
parrot                  Parrotfish data
phyext-methods          A function which returns a phylo4d_ext object
phyextPlot              Plot phylo4d_ext object
phylo4d_ext-class       Class "phylo4d_ext"
plot.censored           Plot a censored tree
plot.taxaset            Plot a taxaset
rcpp_hello_world        Test function
read.analysis.output    Read brownie output to data.frame
read.brownie.string     Read brownie command string into a list
read.characters2        Read CHARACTERS2 block
read.continuous.output
                        Read brownie output from Noncensored analysis
read.discrete.output    Read brownie output from discrete analysis
read.nexus.block        Find a block in a nexus file
read.nexus.simmap       Read a nexus file with SIMMAP-formatted trees
                        (v1.0)
read.ratetest.output    Read brownie output from Censored ratetest
                        analysis
read.simmap             Read simmap-formatted newick strings (version
                        1.0)
read.simmap.new         Read simmap-formatted newick strings (version
                        1.5)
readBrownie             Read in a brownie formatted nexus file
rm.ending               Remove semi-colon from end of a string
run.analysis            Execute a brownie-formatted nexus file
run.asis                Execute a brownie object
runCensored             Add and Run Censored ratetest command to
                        brownie object
runDiscrete             Add and Run a discrete command to brownie
                        object
runNonCensored          Add and Run NonCensored ratetest command to
                        brownie object
scan.textout            Scan raw brownie analysis text
scan.treesout           Scan raw brownie tree text
showBrownie             Show brownie command functions
showSubNodes            Display subnodes
strip                   Strip whitespace
summaryCont             Summarize the results of a Non-censored rate
                        test.
summaryDiscrete         Summarize the results of a discrete character
                        evolution analysis.
summaryRatetest         Summarize the results of a Censored rate test.
taxData                 Return datatype string representing the taxaset
                        datatype
taxa.charvect           Return a character vector of taxa
taxa.mono               Are taxa monophyletic
taxaname.to.taxind      Taxaset index to taxaset name
taxaset.names           Get names of the taxasets
taxaset.rename          Format a taxaset name for use internally
taxind.to.dataind       Convert taxaset column index to data slot
                        column index
validPhylo4d_ext        Validate a phylo4d_ext object
which.mono              Find trees in a list which contain a monophyly
                        for a certain taxaset
which.para              Find trees in a list which contain a paraphyly
                        for a certain taxaset
write.brownie.matrix    Return a matrix formatted for the brownie core
write.brownie.string    Convert commands list to commands string
write.brownie.vector    Return a vector formatted for the brownie core
write.nexus.simmap      Write phylo4d_ext object or list to a nexus
                        file
write.simmap            Write phylo4d_ext tree to newick string
}

Further information is available in the following vignettes:
\tabular{ll}{
\code{getting_started} \tab Getting started with RBrownie (source, pdf)\cr
\code{users_guide} \tab RBrownie User's Guide (source, pdf)\cr
}

}
\author{
J. Conrad Stack <stack@psu.edu>, Brian O'Meara, Luke Harmon

Maintainer: J. Conrad Stack <stack@psu.edu>

}
\references{

}
\keyword{ package }

