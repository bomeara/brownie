#---------------------------------------------
# Methods for manipulating brownie objects
# -	Extra methods for extracting different information
# 	from nexus-formatted files.
#
#---------------------------------------------

# set generics
setGeneric("commands", function(x) { standardGeneric("commands")} )
setGeneric("commands<-", function(x,value) { standardGeneric("commands<-")} )
setGeneric("removeCommands", function(x,index) { standardGeneric("removeCommands")} )
setGeneric("weight", function(x) { standardGeneric("weight")} )
setGeneric("weight<-", function(x,value) { standardGeneric("weight<-")} )
setGeneric("datatypes", function(x) { standardGeneric("datatypes")} )
setGeneric("datatypes<-", function(x,enforce=TRUE,value) { standardGeneric("datatypes<-")} )
setGeneric("taxasets", function(x) { standardGeneric("taxasets")} )
setGeneric("taxasets<-", function(x,value) { standardGeneric("taxasets<-")} )


## OVERLOADED 'phylo4d' functions:
# Risky ....
# overload addData (addData the normal way, then guess it's datatype)
setMethod("addData", signature(x="brownie"),
	function(x,...,known.types=NULL) {
		
		ncols.prev = ncol(tdata(x))
		dtypes.prev = datatypes(x)
		
		# call original phylobase function:
		x = getMethod("addData","phylo4d")(x,...)
		
		# now add datatypes if necessary:
		ncols.now = ncol(tdata(x))
		if(ncols.prev==ncols.now)
			return(x)
		
		newcols = seq(ncols.prev+1,ncols.now)
		
		if(missing(known.types) || is.null(known.types)){
			known.types = .guess.datatype(tdata(x)[,newcols,drop=F])
		} else {
			if(length(known.types) != length(newcols))
				warning("known.types does not match the number of new fields")
		}
		datatypes(x) <- c(dtypes.prev,known.types)
		return(x)
		
})

# Show brownie object:
showBrownie <- function()
{
	cat("brownie class help:\n")
	cat("-------------------------------\n")
	cat("Show/Set datatypes: datatypes(brownieobj) <- c(contData(),discData(),taxaData())\n")
	cat("Show/Set brownie commands: commands(brownieobj) <- c('cmd1','cmd2',....)\n")
	cat("Show/Add taxa sets: taxasets(brownieobj) <- c('taxa1','taxa2',....) \n")
	cat("Show/Add data columns: sndata(brownieobj,<column_index>) <- data.frame\n")
	cat("Remove data columns: rmdata(brownieobj,<column_index>)\n")
	cat("Show which branches subnodes are on: snbranch(brownieobj)\n")
	cat("Show where on branch subnodes are: snposition(brownieobj)\n")
	cat("-------------------------------\n")	
}

setMethod("show","brownie", function(object){ printphylo4(object); showSubNodes(object); showBrownie()})


#---------------
## COMMANDS:
#
#---------------

setMethod("commands", signature(x="brownie"),
  function(x) {
	return(x@commands)
})

# actually this appends... (not a great idea)
setReplaceMethod("commands", signature(x="brownie"),
	function(x,value) {
		
		cmdtext = character(0)
		if(is.list(value)){
			cmdtext = write.brownie.string(value)
		} else {
			cmdtext = value[1]
		}
		
		x@commands = append(x@commands,cmdtext)
		x		
})


hasCommands <- function(x)
{
	if(!is(x,'brownie'))
		return(FALSE)
	
	return((length(commands(x))!=0))
}


setMethod("removeCommands", signature(x="brownie",index="numeric"),
	function(x,index){	
		if(all(length(x@commands) >= index))
		{
			x@commands = x@commands[-index]
		}
		x
})


# index here is the command you want to remove.  All instances of that 
# brownie command are removed
#
setMethod("removeCommands", signature(x="brownie",index="character"),
	function(x,index){
		
		cmdnames = unname(sapply(commands(shit), function(i) read.brownie.string(i)$command ))
		indices = which(tolower(cmdnames) == tolower(index))
		if(length(indices)>0)
			return(removeCommands(x,indices))
			
		return(x)
})



#---------------
## WEIGHT:
#
#---------------
setMethod("weight", signature(x="brownie"),
  function(x) {
	return(x@weight)
})

setReplaceMethod("weight", signature(x="brownie"),
  function(x,value) {
	x@weight = value
	return(x)
})




#---------------
## DataTypes
#
#---------------

setMethod("datatypes", signature(x="brownie"),
  function(x) {
	return(x@datatypes)
})


setReplaceMethod("datatypes", signature(x="brownie"),
  function(x,enforce=TRUE,value) {
	if(length(value) == ncol(tdata(x)) || !enforce){
		
		if(enforce)
		{
			if(all(value %in% brownie.datatypes())){
				x@datatypes = value
			} else {
				warning("Some values are not in the regular set (",paste(brownie.datatypes(),collapse=","),")")
			}
		}else{
			x@datatypes = as.character(value)
		}
	
	} else {
		warning("Values was not the same length as the number of data columns")
	}
	return(x)
})


setMethod("rmdata", signature(x="brownie",index="numeric"),
  function(x,index) {
	  getMethod("rmdata",signature("phylo4d_ext","numeric"))(x,index)
		if(length(index)>0 && index <= ncol(tdata(x)))
		{
			x@data = x@data[,-index,drop=F]
			dtypes = datatypes(x)
			datatypes(x) <- dtypes[-index]
		}
	return(x)
})




#---------------
## TAXASETS:
#
#---------------

setMethod("taxasets", signature(x="phylo4d"),
  function(x) {
	retdat = data.frame(NULL)
	dat = tdata(x,"tip")
	colinds = grep("^TAXSET_",names(dat))
	if(length(colinds) > 0)
		retdat = dat[,colinds,drop=F]
	return(retdat)
})


setMethod("taxasets", signature(x="brownie"),
  function(x) {
	
	retdat = data.frame(NULL)
	dat = tdata(x,"tip")
	
	colinds = which(datatypes(x) == taxaData())
	if(length(colinds) > 0)
		retdat = dat[,colinds,drop=F]
	
	return(retdat)
})


hasTaxasets<-function(x)
{
	retbool = (ncol(taxasets(x))!=0)
	retbool
}

# adds only:
setReplaceMethod("taxasets", signature(x="phylo4d"),
  function(x,value) {
	
	setname = FALSE
	
	# if it is a vector and not a
	if(is.null(dim(value)))
	{
		
		if(is.character(value))
		{
			# assume it is a character vector of the new taxa to add
			ff = factor(value,levels=tipLabels(x))
			if(all(table(ff)==0))
				stop("Character vector provided does not contain any valid taxa names")
			
			value = data.frame(as.character(table(ff)),row.names=as.character(levels(ff)))
			
		} else {
			# assume binary data with names which are the taxa
			value = data.frame(value)
		}
		colnames(value) <- "TAXSET_"
	}
	
	if( !all(sort(rownames(value)) == sort(tipLabels(x))) )
	{
		stop("All taxa in assigning value need to equate to taxa present")
	} else {
		
		if(is.null(colnames(value)))
			colnames(value) <- "TAXSET_"
		
		# ohh, this is a hack:
		if(is(x,"brownie")){
			x <- addData(x,value,known.types=taxaData())
		} else {
			x <- addData(x,value)
		}
	}
	
	return(x)
})


setReplaceMethod("taxasets", signature(x="brownie"),
  function(x,value) {
	if(!existsMethod("taxasets<-","phylo4d"))
		stop("Taxaset replacement method not available: package bug, email authors")
	
	x = getMethod("taxasets<-","phylo4d")(x,value)
	#datatypes(x) <- append(datatypes(x),taxaData()) # comments out after hack applied
	
	return(x)
})


# Internal:
# Return character vector of taxa from taxasets
#
taxa.charvect <- function(x,taxind)
{
	if(!hasTipData(x))
		stop("x does not have any data")

	retbool = FALSE
	tset = NULL
	tsets = taxasets(x)
	if(is.character(taxind[1]))
	{
		# assume it's a column header
		tset = tsets[taxind]
	} else {
		tset = tsets[,taxind[1],drop=F]
	}
	tipvect = rownames(tset)[which(tset == 1)]

	return(tipvect)
}

# Internal:
# Is a character vector of taxa monophyletic on phylogeny x?
#
taxa.mono <- function(x,tipvector)
{
	retbool = FALSE
	desc = descendants(x,MRCA(x,tipvector))
	if(length(desc) == length(tipvector))
		retbool = all(sort(names(desc)) == sort(tipvector)) # this should always be true I think
	
	retbool
}

# check if certain taxaset is monophyletic:
#
areTaxaMono <- function(x,taxindex)
{	
	return(taxa.mono(x,taxa.charvect(x,taxindex)))
}

# Check if taxaset is paraphyletic
#
areTaxaPara <- function(x,taxind)
{	
	retbool = FALSE
	tipvect = taxa.charvect(x,taxind) # character vector of taxa to check
	
	desc = descendants(x,MRCA(x,tipvect))
	excluded = setdiff(names(desc),tipvect)
	if(length(excluded) > 0)
		retbool = taxa.mono(x,excluded)
	
	return (retbool)
}

