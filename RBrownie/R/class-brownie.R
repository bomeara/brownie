#--------------------------------------
# Brownie class
# @author Conrad Stack
#
# - This class wraps a phylo4d_ext object.
# 	It also includes brownie commands,
# 	which are used to run the analysis.
#--------------------------------------


# @param phylo4d a phylogeny with data
# @param commands text commands to be placed in a brownie block
# @param datatypes annotation for the data part of phylo4d
setClass("brownie",
		representation("phylo4d_ext",commands="character",datatypes="character")
		)


#---------------
# constructors |
#---------------

# 'brownie' to be a generic function
setGeneric("brownie", function(x, ...) { standardGeneric("brownie")} )


## first arg is a phylo
setMethod("brownie", "phylo", function(x, annote=list()){
	
	converted = as(x,"brownie")
	converted@annote <- annote
	
	return(converted)
})



#-------------------------------
# Validate correct brownie object
#-------------------------------

# use phylobase checker
checkBrownie <- function(object)
{
	retval = TRUE
	
	# check brownie commands if needed
	if(length(object@commands)!=0 && !all(object@commands==""))
		retval = retval && checkCommands(object@commands)
	
	if(length(object@datatypes)!=0 && !all(object@datatypes==""))
		retval = retval && checkDataTypes(object@datatypes)
	
	return (retval && checkPhylo4(object))
}
setValidity("brownie",checkBrownie)


# TODO: fill in with stuff from processBrownie.R
checkCommands <- function(cmds)
{
	return (TRUE)
}

checkDataTypes <- function(dtypes)
{
	
	# These need to be either: "taxset","cont","discrete","data"
	# brownie_datatype_options should be "more" global
	# 
	return (all(dtypes %in% brownie.datatypes()))
}

