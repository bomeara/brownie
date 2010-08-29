\name{addDiscrete}
\alias{addDiscrete}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add Discrete character models and reconstructions to brownie object
}
\description{
This function adds a character string to the the commands slot of a brownie object (or list of brownie objects) which instructs the brownie program to run a discrete character reconstruction.
}
\usage{
addDiscrete(obj, 
	file = NULL, 
	model = brownie.models.discrete()[1], 
	model.state = NULL, 
	ratemat = NULL, 
	freq = brownie.freqs()[1], 
	statevector = NULL, 
	treeloop = FALSE, 
	charloop = FALSE, 
	allchar = FALSE, 
	variable = FALSE, 
	reconstruct = FALSE, 
	breaknum = 0, 
	fileappend = FALSE, 
	filereplace = FALSE, 
	globalstates = FALSE, 
	usetempfile = FALSE, 
	...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
A brownie object (or list of brownie objects)
}
  \item{file}{
The brownie core write the output of the analysis to this file (optional)
}
  \item{model}{
Allows you to specify whether to use a USER-specified model, a model where all rates are EQUAL, a REVersible model where q_ij=q_ji for all states i and j but are otherwise free to vary, or an NONREVersible model where all rates can vary independently. (see brownie.models.discrete() for list of options)
}
  \item{model.state}{
Which model state would you like to use?
}
  \item{ratemat}{
If model == USER, then this is a vector containing information about rate parameters.  The model specification, except for the built-in types (equal, rev, nonrev), is grossly similar to PAUP's method for specifying which rates are constrained to be equal, plus allows fixing of certain values. For example, the following rate matrix:

                to->  0     1     2
               from ------------------
                 0 |  -     a    0.5
                 1 |  b     -     c
                 2 | 0.0    a     -
                 
means that rate q01 (instantaneous rate going from state 0 to state 1), with value a, must also equal rate q21 but is otherwise unconstrained (so they are both optimized, but forced to take the same value), while q20 is forced to a rate value of 0.0. We could have specified any non-negative fixed value. Basically, all rates sharing a letter take the same optimized rate value, while those assigned a number have a fixed rate value. To specify the above rate matrix, the command would be "discrete model=user ratemat=(a 0.5 b c 0.0 a)"
(see write.brownie.matrix for help formatting)
}
  \item{freq}{
The probability of each state at the root can be based on the EMPIRICAL distribution at the tips, can be SET by the user (using the statevector command), can be OPTIMIZEd as part of the model, can be set to EQUILIBRIUM frequencies (the frequencies expected with the optimized rate matrix given infinitely-long branches), or can be set to be UNIFORM (equal).
}
  \item{statevector}{
Contains the user-specified probabilities of the ancestral states. Example: "discrete freq=user statevector=(0.4 0.6)" for a binary trait.  
(see write.brownie.vector for help formatting)
}
  \item{treeloop}{
Should the analysis loop across the trees?  To use a specific tree then run the addChoose(tree=<TREENUM>) before this command.
}
  \item{charloop}{
Should the analysis loop across characters?  To use a specific character then run the addChoose(char=<CHARNUM>) before this command.
}
  \item{allchar}{
Estimates the model parameters using all the characters simultaneously, not based on a single character.
}
  \item{variable}{
If true, performs a correction to correct rates for only examining variable characters.
}
  \item{reconstruct}{
Using the likelihood rate matrix and state frequencies, reconstructs the joint estimates of the ancestral states at internal nodes and, optionally, along the branches. This uses the Pupko et al 2000 algorithm, which is fast but only returns the estimated states, not the confidence in these states.
}
  \item{breaknum}{
Setting this value >0 allows the program to estimate the likeliest state at breaknum points along each branch. This can be useful in estimating when on a branch a character changed. Note that this may underestimate the number of changes on a branch (for example, on a very long branch which starts and ends in state 0, with high enough transition rates it may be probable that the character has changed  multiple times on that branch, but at any given instant on that branch, the likeliest state it will be in is state 0).
}
  \item{fileappend}{
If the output file exists, appends to it rather than overwrites it. Turned on by default.
}
  \item{filereplace}{
If set to yes, if the output file already exists it will be quietly replaced.
}
  \item{globalstates}{
If no, the number of character states assumed for each character is the maximum number of observed states for just that one character. If yes, the number of states for each character is the maximum number of states observed for any character, even if the observed character is lacking some of those states. This is useful if, for example, you have simulated a three state character on the tree for parametric bootstrapping but some of the simulations result in characters with just states 0 and 1. Globalstates will automatically be set to yes if allchar=y
}
  \item{usetempfile}{
Should a tempfile be written to (instead of file).  Not currently used.
}
  \item{\dots}{
Not currently used
}
}
\details{
This is a function to calculate the likelihood estimates of discrete character evolution parameters and scores for these models. This will allow you to do things like compare models with equal gain and loss rates with models which allow these to vary, evaluate models with a mixture of fixed and free rates, reconstruct the joint likelihood estimates of ancestral states at nodes and at various points within branches, and more.

From Brownie command-line manual:
Usage: Discrete [model=] [freq=] [ratemat=] [statevector=] [treeloop=] [charloop=] [allchar=] [variable=] [reconstruct=] [breaknum=] [file=] [append=] [replace=] [globalstates=]
}
\value{
The brownie object \code{obj} with a new command added.
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
data(geospiza_ext)
brownierun=runDiscrete(geospiza_ext,brfile="disc_junk.txt",
models=c("nonrev","nonrev"),freqs=c("unif","equilib"),reconstruct=TRUE)
brownierun
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

