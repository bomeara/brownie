
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

# TODO
checkCommands <- function(cmds)
{
	return (TRUE)
}

checkDataTypes <- function(dtypes)
{
	
	# These need to be either: "taxset","cont","discrete","data"
	# brownie_datatype_options should be "more" global
	# 
	return (all(dtypes %in% brownie_datatype_options))
	
}

