#--------------------------------------------------
# Supplemental brownie functions:
# -	Need to load these things first
#
#--------------------------------------------------

## Valid 'Datatypes':
#
# for checking datatypes (to be used with 'data' slot)
brownie.datatypes <-function()
{
	brownie_datatype_options = c("taxset","cont","discrete","undef")
	return(brownie_datatype_options)
}

discData <- function() { return("discrete") }
contData <- function() { return("cont") }
taxaData <- function() { return("taxset") }
genericData <- function() { return("undef") }



## Valid Brownie Commands:
#
#
brownie.commands <- function(with.desc=F)
{
	
	brownie_cmd_options = c(
		"cont",
		"ratetest",
		"discrete",
		"model",
		"choose"
	)
	
	brownie_cmd_description = c(
		"Non-censored rate test",
		"Censored rate test",
		"Discrete ancestral state reconstruction",
		"Specifying model type",
		"Specify with (tree|...) to use"
	)
	
	if(with.desc)
		brownie_cmd_options = cbind(brownie_cmd_options,brownie_cmd_description)
	
	return(brownie_cmd_options)
}


## Valid brownie models
#
brownie.models <- function(with.desc=F)
{
	brownie_model_options = c(
		"BM1",
		"BMS",
		"BMC",
		"OUSM",
		"OUCM"
	)
	
	brownie_model_description = c(
		"Brownian motion, one rate parameter",
		"Brownian motion, with different rate parameters for each state on a tree",
		"Brownian motion, with one rate parameter for branches with state changes and another for branches without change.",
		"Ornstein-Uhlenbeck with one mean per discrete state (attraction and rate parameters constant across tree)",
		"Ornstein-Uhlenbeck with independent means on branches with and without changes in a discrete character (attraction and rate parameters constant across tree)"
	)
	
	if(with.desc)
		brownie_model_options = cbind(brownie_model_options,brownie_model_description)
	
	return(brownie_model_options)
}


## Valid Brownie options (used with cmds)
#
#	Note: 	this will usually need to be called twice; once to 
#			get the option names and another to get the validate
#			functions:
#			
brownie.options <- function(with.decs=F, with.vals=F)
{
	
	brownie_option_options = c(
		"taxset",
		"treeloop",
		"charloop",
		"reps",
		"tree",
		"type",
		"file",
		"append",
		"replace"
	)
	
	brownie_option_description = c(
		rep("",brownie_option_options)
	)
	
	brownie_option_validate = c(
		checkval.dummy,
		checkval.yesno,
		checkval.yesno,
		checkval.integer,
		checkval.integer,
		checkval.model,
		checkval.dummy,
		checkval.yesno,
		checkval.yesno
	)
	
	if(with.desc)
		brownie_option_options = cbind(brownie_option_options,brownie_option_description)
	
	# if this, then only return the function list
	if(with.vals){
		names(brownie_option_validate)<-brownie_option_options
		brownie_option_options<-brownie_option_validate
	}
	
	return(brownie_option_options)
}


# check option values:
checkval.dummy <-function(optstr)
{
	return(TRUE)
}

checkval.yesno <-function(optstr)
{
	return( as.logical(length(grep("(y|yes|n|no)",tolower(optstr)))) )
}

checkval.integer <- function(optstr)
{
	return( is.na(as.integer(optstr)) )
}

checkval.model <- function(optstr)
{
	# TODO: figure out if brownie is in fact case-insensitive like I'm assuming
	#
	return( tolower(optstr) %in% tolower(brownie.models) )
}



