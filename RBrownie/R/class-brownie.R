setClass("brownie",
		representation("phylo4d",commands="character",datatypes="character"),
		validity = checkBrownie
		)

		
#---------------
# constructors |
#---------------

# for checking datatypes (to be used with 'data' slot)
brownie_datatype_options = c("taxset","cont","discrete","data")

# 'brownie' to be a generic function
setGeneric("brownie", function(x, ...) { standardGeneric("brownie")} )


## first arg is a phylo
setMethod("brownie", "phylo", function(x, annote=list()){
	
	converted = as(x,"brownie")
	converted@annote <- annote
	
	return(converted)
})



