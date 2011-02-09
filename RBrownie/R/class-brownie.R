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
		representation(commands="character",datatypes="character"),
		prototype=prototype(commands="",datatypes=""),
		contains="phylo4d_ext"
		)


#---------------
# constructors |
#---------------

# 'brownie' to be a generic function
setGeneric("brownie", function(x, ...) { standardGeneric("brownie")}, valueClass=c("brownie","list") )


## first arg is a phylo
setMethod("brownie", "phylo", function(x, annote=list()){
	
	# special way to convert phylo to brownie
	converted = as(x,"brownie")
	converted@annote <- annote
	tipLabels(converted) <- checkLabel(tipLabels(converted))
	
	return(converted)
})

setMethod("brownie", "phylo4", function(x, annote=list()){
	
	converted = as(x,"phylo4d")
	return(brownie(converted))
})

setMethod("brownie","phylo4d", function(x,dataTypes,annote=list()){
	
	if(hasSingle(x))
		x <- phyext(x)
	
	if(missing(dataTypes)){
		converted = new("brownie",x,datatypes=.guess.datatype(tdata(x)))
	} else {
		converted = new("brownie",x,datatypes=dataTypes)
	}
	
	# TODO: -check if datatypes match up with the @data slot, 
	#		 convert datatypes if necessary.

	# Convert to valid names (if necessary)
	tipLabels(converted) <- checkLabel(tipLabels(converted))	
	
	return(converted)
})

setMethod("brownie","list",
	function(x,...){
		x = sapply(x,brownie,...)
		return(x)
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


# Internal function to parse / check assumptions block
.process.assumptions <- function(obj,block.txt,simbool)
{
	
	if(!is(obj[[1]],"brownie"))
		stop("Processed object needs to be of class brownie")
	
	for(aline in block.txt)
	{
		tokens = strsplit(aline,"\\s")[[1]]
		tokens = .split.tokens(tokens,"=")
		next.is.name=F
		next.is.taxa=F
		taxinds=numeric(0)
		taxname=""
		
		for(kk in seq(length(tokens)))
		{
			if(tolower(tokens[kk]) == "taxset")
			{
				next.is.name = T
				next
			} else {
				
				# reset and start recording names
				if(next.is.name){
					next.is.name=F
					taxname = sub("=","",tokens[kk])
					next.is.taxa=T
					next
				}
				
				if(next.is.taxa){
					if(length(grep(";",tokens[kk]))==1)
					{
						taxinds = append(taxinds,sub(";","",tokens[kk]))
						break
					} else {
						taxinds = append(taxinds,tokens[kk])
					}
				}
			}
		}
		
		# add this taxaset to the object, if there is one:
		if(taxname != "")
		{
			taxinds.as.range = FALSE
			nm = paste("TAXSET",taxname,sep="_")
			taxaI = data.frame(all=rep(0,nTips(obj[[1]])))
			names(taxaI) <- nm
			
			# check for taxinds as numbers:
			if(length(taxinds) == 1)
			{
				both = get.left.right(taxinds,"-")
				if( length(both) == 2 && !any(is.integer(both)) ){
					taxinds.as.range = TRUE
					taxinds = as.integer(both)
				}
			}
			
			for(tind in seq(length(obj)))
			{
				# convert if needed
				if(!inherits(obj[[tind]],"phylo4d"))
					obj[[tind]] = phylo4d(obj[[tind]])
				
				if(!taxinds.as.range){
					
					if(length(grep("_",taxinds))!=0)
						warnings("Stipping underscores (because phylobase does).")
					
					# assume that phylobase will remove underscore characters
					#if(!simbool)
					#	taxinds = gsub("_","",taxinds)
					
					taxaI[,1] = sapply(tipLabels(obj[[tind]]),function(i) ifelse(i %in% taxinds,1,0),simplify=T)
				} else {
					taxaI[seq(taxinds[1],taxinds[2]),1] = taxaI[seq(taxinds[1],taxinds[2]),1] = 1
				}
				#names(taxaI) <- nm
				obj[[tind]] = addData(obj[[tind]],tip.data=taxaI)
				
			}
		}
	}
	
	return (obj)
}


# Internal:
is.binary <- function(seqvec)
{
	
}

# TODO: Make sure is.factor should be used here to select discrete data
# Internal:
.guess.datatype <- function(datvals)
{
	if(is.null(dim(datvals)))
		return(character(0))
	
	ndatcols = ncol(datvals)
	datatypes = rep(genericData(),ndatcols)
	if(ndatcols > 0)
	{
		datatypes[sapply(seq(ndatcols),function(i) is.factor(datvals[,i]))] = discData()
		datatypes[sapply(seq(ndatcols),function(i) is.numeric(datvals[,i]))] = contData()
		datatypes[grep("TAXSET_",names(datvals))] = taxData()
	}
	
	return(datatypes)
}


# Internal function to parse / check assumptions block
.process.datatypes <- function(obj)
{
	if(!is(obj[[1]],"brownie"))
		stop("Processed object needs to be of class brownie")
	
	datvals = tdata(obj[[1]])
	ndatcols = ncol(datvals)
	if(ndatcols > 0)
	{
		dtypes = .guess.datatype(datvals)
		
		for(tind in seq(length(obj)))
		{
			datatypes(obj[[tind]]) <- dtypes
		}
	}
	
	return(obj)
}



# internal function to parse / check brownie block
.process.brownie <- function(obj,block.txt)
{
	# TODO
	return(obj)
}

