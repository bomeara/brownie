#---------------------------------------------
# Methods for manipulating brownie objects
# -	Extra methods for extracting different information
# 	from nexus-formatted files.
#
#---------------------------------------------

# set generics
setGeneric("commands", function(x) { standardGeneric("commands")} )
setGeneric("commands<-", function(x,add=FALSE,index=NULL,replace=FALSE,value) { standardGeneric("commands<-")} )
setGeneric("clearCommands", function(x) { standardGeneric("clearCommands")} )
setGeneric("removeCommands", function(x,index) { standardGeneric("removeCommands")} )
setGeneric("datatypes", function(x) { standardGeneric("datatypes")} )
setGeneric("datatypes<-", function(x,enforce=TRUE,value) { standardGeneric("datatypes<-")} )
setGeneric("taxasets", function(x) { standardGeneric("taxasets")} )
setGeneric("taxasets<-", function(x,taxnames,value) { standardGeneric("taxasets<-")} )
setGeneric("removeTaxasets", function(x,index) { standardGeneric("removeTaxasets")} )
setGeneric("hasTaxasets", function(x) { standardGeneric("hasTaxasets")} )






## OVERLOADED 'phylo4d' functions:
#
setMethod("addData",signature(x="list"),
	function(x,...){
		x = sapply(x,addData,...)
	return(x)
})

# overload addData (addData the normal way, then guess it's datatype)
setMethod("addData", signature(x="brownie"),
	function(x,...,dataTypes=NULL) {
		
		ncols.prev = ncol(tdata(x))
		dtypes.prev = datatypes(x)
		
		# call original phylobase function:
		x = getMethod("addData","phylo4d")(x,...)
		
		# now add datatypes if necessary:
		ncols.now = ncol(tdata(x))
		if(ncols.prev==ncols.now)
			return(x)
		
		newcols = seq(ncols.prev+1,ncols.now)
		
		if(missing(dataTypes) || is.null(dataTypes)){
			dataTypes = .guess.datatype(tdata(x)[,newcols,drop=F])
		} else {
			if(length(dataTypes) != length(newcols))
				warning("dataTypes does not match the number of new fields")
		
			
			# Added (7/26/2010)
			# Making sure that datatypes is either a factor or numeric vector
			#
			for(ii in seq(length(newcols)))
			{
				dat = tdata(x)[,newcols[ii],drop=T]
				
				if(dataTypes[ii] == contData() && !is.contData(dat))
					tdata(x)[,newcols[ii]] <- as.contData(dat)
				
				if(dataTypes[ii] == discData() && !is.discData(dat))
					tdata(x)[,newcols[ii]] <- as.discData(dat)
				
			}
			
		}
		datatypes(x) <- c(dtypes.prev,dataTypes)
		return(x)
})

# Show brownie object:
showBrownie <- function()
{
	cat("brownie class help:\n")
	cat("-------------------------------\n")
	cat("Show/Set datatypes: datatypes(brownieobj) <- c(contData(),discData(),taxData())\n")
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

setMethod("commands", signature(x="list"),
  function(x) {
	return(x[[1]]@commands)
})

# actually this appends... (not a great idea)
setReplaceMethod("commands", signature(x="brownie"),
	function(x,add=FALSE,index=NULL,replace=FALSE,value) {
		
		cmdtext = character(0)
		if(is.list(value)){
			cmdtext = write.brownie.string(value)  # assume that list is a 'command struct'-type list
		} else {
			cmdtext = value
		}
		

		if(add)
		{
			if(is.null(index))
			{
				x@commands = append(x@commands,cmdtext)
			} else {
				if(max(index) > length(x@commands))
					stop("Index is too large: cmdlen = ",length(x@commands),"; indices = ",index)
				
				if(length(index) != length(cmdtext))
					stop("replacement index and value do not have the same length")
				
				tmp = x@commands[index]
				x@commands[index] = value
				if(!replace)
				{
					if(max(index) == length(x@commands)){
						x@commands = append(x@commands,tmp)
					} else {
						tmp2 = x@commands[seq(from=max(index)+1,to=length(x@commands))]
						x@commands = c(x@commands[1:max(index)],tmp,tmp2)
					}
				}
				
			}
		} else {
			x@commands = cmdtext
		}
		x		
})


setMethod("clearCommands",signature(x="brownie"),
	function(x){
		x@commands <- character(0)
	return(x)
})


setMethod("clearCommands",signature(x="list"),
	function(x){
		x = sapply(x,clearCommands)
	return(x)
})


hasCommands <- function(x)
{
	if(!is.list(x))
		x = list(x)
	
	if(!is(x[[1]],'brownie'))
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
		
		cmdnames = unname(sapply(commands(x), function(i) read.brownie.string(i)$command ))
		indices = which(tolower(cmdnames) == tolower(index))
		if(length(indices)>0)
			return(removeCommands(x,indices))
			
		return(x)
})


setMethod("removeCommands",signature(x="list",index="ANY"),
	function(x,index){
		
		x = sapply(x,removeCommands,index)
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


setMethod("datatypes", signature(x="list"),
  function(x) {
	return(x[[1]]@datatypes)
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


setReplaceMethod("datatypes", signature(x="list"),
  function(x,enforce=TRUE,value) {
	  for(ii in seq(length(x)))
	  	datatypes(x[[ii]],enforce=enforce) <- value
	  
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
	
	colinds = which(datatypes(x) == taxData())
	if(length(colinds) > 0)
		retdat = dat[,colinds,drop=F]
	
	return(retdat)
})

setMethod("taxasets",signature(x="list"),
	function(x) {
		return(taxasets(x[[1]]))
})

setMethod("hasTaxasets",signature(x="phylo4d"),
	function(x)
	{
		retbool = (ncol(taxasets(x))!=0)
		retbool
})

setMethod("hasTaxasets",signature(x="ANY"),
	function(x)
	{
		return(FALSE)
})

setMethod("hasTaxasets",signature(x="list"),
	function(x){
		return(hasTaxasets(x[[1]]))
})



# adds only:
setReplaceMethod("taxasets", signature(x="phylo4d"),
  function(x,taxnames,value) {
		
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
			value = data.frame(value,row.names=tipLabels(x))
		}
		
		if(missing(taxnames)){
			colnames(value) <- "TAXSET_"
		}
	}
	
	if( !all(sort(rownames(value)) == sort(tipLabels(x))) )
	{
		stop("All taxa in assigning value need to equate to taxa present")
	} else {

		if(!missing(taxnames))
			colnames(value) <- taxnames
					
		cnames = colnames(value)
		notlabeled = grep("^[^T][^A][^X][^A][^S][^E][^T][^_]",cnames)
		if(length(notlabeled) > 0){
			cnames[notlabeled] = paste("TAXSET_",cnames[notlabeled],sep="")
			colnames(value) <- cnames
		}

		# ohh, this is a hack:
		if(is(x,"brownie")){
			x <- addData(x,value,dataTypes=taxData())
		} else {
			x <- addData(x,value)
		}
	}
	
	return(x)
})


setReplaceMethod("taxasets", signature(x="brownie"),
  function(x,taxnames,value) {
	if(!existsMethod("taxasets<-","phylo4d"))
		stop("Taxaset replacement method not available: package bug, email authors")
	
	x = getMethod("taxasets<-","phylo4d")(x,taxnames,value)
	#datatypes(x) <- append(datatypes(x),taxData()) # comments out after hack applied
	
	return(x)
})

setReplaceMethod("taxasets",signature(x="list"),
	function(x,taxnames,value) {
	for(ii in seq(length(x)))
		taxasets(x[[ii]],taxnames=taxnames) <- value
	
	return(x)
})


setMethod("removeTaxasets", signature(x="brownie",index="character"), 
	function(x,index) {
		tnames = names(taxasets(x))
		altnames = sub("^TAXSET_(.*)$","\\1",tnames)
		if(index %in% tnames)
		{
			return(removeTaxasets(x,index = which(index == tnames)[1]))
		} else {
			if(index %in% altnames)
			{
				return(removeTaxasets(x,index = which(index == altnames)[1]))
			} else{
				stop("Could not find taxaset with name ",index)
			}
		}
})


setMethod("removeTaxasets", signature(x="brownie",index="numeric"), 
	function(x,index) {
	datind = taxind.to.dataind(br,index)
	if(length(datind) != 1)
		stop("problem converting taxa index into data index")
	
	br@data = br@data[,-datind,drop=F]
	br@datatypes = br@datatypes[-datind]
	return(br)
})

	

#
taxaname.to.taxind <- function(x,taxnames)
{
	taxinds = integer(length(taxnames))
	for(ii in seq(length(taxnames)))
	{
		taxinds[ii] = which(tipLabels(x) == taxnames[ii])
	}
	return(taxinds)
}

#
taxind.to.dataind <- function(x,taxind)
{
	index = integer(0)
	if(is.character(taxind)){
		index = which(colnames(tdata(x,'tip'))==taxind)
		if(length(index) == 0)
			index = which(colnames(tdata(x,'tip'))==taxaset.rename(taxind))
		
		if(length(index) == 0)
			stop("Could not find taxaset called: ",taxind,"\n These are available:",taxaset.names(x))
	} else {
		if(taxind <= length(taxasets(x))){
			index = which(colnames(tdata(x,'tip')) == colnames(taxasets(x))[taxind])
		}else{
			stop("Index ",taxind," is out of range.")
		}
	}
	return(index)
}


# get taxaset names
taxaset.names <- function(x)
{
	retnames = character(0)
	
	if(hasTaxasets(x))
	{
		retnames = colnames(taxasets(x))
		retnames = sub("TAXSET_","",retnames)
	}
	return(retnames)
}

# append TAXSET_ to front of x
taxaset.rename <- function(x)
{
	return(sprintf("TAXSET_%s",x))
}

# check if taxaset has "TAXSET_" prefix:
is.taxaname.internal <- function(x)
{
	return( (length(grep("TAXSET_",x)) > 0) )
}

# Internal:
# Return character vector of taxa from taxasets
# NOTE: -If taxind is a character string and does not have an 'internal representation' 
#		 (prefixed with TAXSET_), then TAXSET_ is added.
#
taxa.charvect <- function(x,taxind,append.internal=TRUE)
{
	if(!hasTipData(x))
		stop("x does not have any data")

	retbool = FALSE
	tset = NULL
	tsets = taxasets(x)
	if(is.character(taxind[1]))
	{
		if(append.internal && !is.taxaname.internal(taxind))
			taxind = taxaset.rename(taxind)
		
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
# with.respect.to is a work in progress; (see ratetest ? for more details)
#
areTaxaPara <- function(x,taxind,with.respect.to=NULL)
{	
	retbool = FALSE
	tipvect = taxa.charvect(x,taxind) # character vector of taxa to check
	tipcompare=character(0)
	if(!is.null(with.respect.to))
		tipcompare = taxa.charvect(x,with.respect.to) 
	
	desc = descendants(x,MRCA(x,tipvect))
	excluded = setdiff(names(desc),tipvect)
	if(length(excluded) > 0)
		retbool = taxa.mono(x,excluded)
	if(!is.null(with.respect.to))
		retbool = retbool && all(sort(excluded) == sort(tipcompare))
	return (retbool)
}

# check if taxasets are mutually exclusive
areTaxaMutex <- function(x,...)
{
	tsets = list(...)
	overlap = character(0)
	
	if(length(tsets) == 0)
		stop("No taxa specified")
	
	if(length(tsets) == 1)
		return(TRUE)
	
	
	for(ii in seq(length(tsets)))
	{
		overlap = append(overlap,taxa.charvect(x,tsets[[ii]]))
	}
	tabo = table(overlap)
	
	return( all(tabo==1) )
	
}

# indices for trees in which the taxaset specified is monophyletic
#
which.mono <- function(treeslist,taxaset)
{
	retints = integer(0)
	if(hasTaxasets(treeslist[[1]]))
	{
		retints = which(sapply(treeslist,areTaxaMono,taxaset,simplify=T))
	} else {
		warning("treeslist does not contain any taxasets")
	}
	return(retints)
}


# Check how many trees in a list are monophyletic for a given taxaset:
checkMono <- function(treeslist,taxaset,percent=T)
{
	if(!is.list(treeslist))
		treeslist = list(treeslist)
	
	indices = which.mono(treeslist,taxaset)
	retnum = length(indices)
	if(percent)
		retnum = retnum / length(treeslist)
	
	return(retnum) 
}



# indices for trees in which the taxaset specified is monophyletic
#
which.para <- function(treeslist,taxaset)
{
	retints = integer(0)
	if(hasTaxasets(treeslist[[1]]))
	{
		retints = which(sapply(treeslist,areTaxaPara,taxaset,simplify=T))
	} else {
		warning("treeslist does not contain any taxasets")
	}
	return(retints)
}


# Check how many trees in a list are paraphyletic for a given taxaset:
checkPara <- function(treeslist,taxaset,percent=T)
{
	if(!is.list(treeslist))
		treeslist = list(treeslist)
	
	indices = which.para(treeslist,taxaset)
	retnum = length(indices)
	if(percent)
		retnum = retnum / length(treeslist)
	
	return(retnum) 
}


## Generic Overloads from phylo4d_ext:

# NOTE: These rm functions assume that datatypes is matched up with tdata(x) and sndata(x) perfectly
setMethod("rmdata", signature(x="brownie",index="numeric",subindex="numeric"),
  function(x,index,subindex) {
	
	if(length(index)>0 && abs(index) <= ncol(tdata(x)) && index != 0)
	{
		x@data = x@data[,-index,drop=F]
		x@datatypes = x@datatypes[-index]
	} else {
		stop("The data column to be removed could not be found in tdata(x)")
	}

	if( hasSubNodes(x) && length(subindex) >0 && abs(subindex) <= ncol(sndata(x)) && subindex != 0)
	{
		x@subnode.data = x@subnode.data[,-index,drop=F]
	}
	
	return(x)
})
	

setMethod("rmdata", signature(x="brownie",index="numeric",subindex="missing"),
  function(x,index) {
	return( rmdata(x,index=index,subindex=index) )
})

setMethod("rmdata", signature(x="brownie",index="character",subindex="missing"),
  function(x,index) {
	ind = which(colnames(tdata(x)) %in% index)
	if(hasSubNodes(x)){
		subind = which(colnames(sndata(x)) %in% index)
	} else {
		subind = numeric(0)
	}
	
	if( !all(sort(ind) == sort(subind)) )
	{
		warning("Problem: datatypes and (@data, @subnode.data) have become out of sync")
	}

	if( length(ind) == 0 && length(subind) == 0 )
	{
		warning("The data column ",index," was not found in brownie object and request to delete it is being ignored")
		return(x)
	} else {
		return( rmdata(x,index=ind,subindex=subind) )
	}
})


