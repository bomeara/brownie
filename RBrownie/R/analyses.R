#
#--------------------------------------------------
# Phylogenetic analyses for RBrownie
# -
#--------------------------------------------------

# run - will run brownie block as-is.  Need to specify what you want returned.
# Ancestral State Reconstruction - 
# Censored Rate tests
# Non-censored Rate tests
# ?




# Method to do ancestral state reconstruction in brownie
#
# return  
#
.run.asr <- function(filename)
{
	# check argument for character string
	if( !is.character(filename) && length(filename) == 1)
		stop(paste("first argument must be a character string specifying one file only"))
	
	# check if filename exists
	if( !file.exists(filename) )
		stop( paste("filename,",filename,", does not exist") )
	
	# convert to unix-style path
	filename = gsub("\\\\","/",filename)
	
	# add single quotes (brownie requirement for full paths
	filename = paste("'",filename,"'",sep="")
	
	return( .Call("doASR",filename,PACKAGE="RBrownie") )
}



# Do ancestral state reconstruction
#
# @param return phylo4d_ext format or non-extended phylo4d
#
asr <- function(phytree,file="",usedata,extended=TRUE)
{
	outtree=NULL
	fname=""
	if(missing(phytree))
	{
		if(file=="")
			stop("Could not locate the file to read")
		fname = file
	} else {
		fname = tempfile()
		if(!write.nexus.both(phytree,fname))
			stop("Writing to temporary file failed")
	}
	
	outtext = .run.asr(fname)[[1]]
	outtext = gsub("'","",outtext)  # brownie seems to put little ticks where the internal node names are
	if(is.simmap(text=outtext)){
		outtree = read.simmap(text=outtext)
		outtree = phyext(outtree)
	} else {
		outtree = read.tree(text=outtext)
	}
	return( outtree )
}

###### Add to commands list:
#
addcmd.taxaset <- function(obj,cmdobj,x)
{
	# set 'taxset' option:
	tsnames = taxaset.names(obj)
	if(is.character(x))
	{
		# TODO: add new taxaset if it doesn't exist.
		if(x != "ALL" && !(x %in% tsnames))
		{
			stop("Specified taxaset,",x,", is not a defined taxaset.")
		} else {
			cmdobj$options = rbind(cmdobj$options,c("taxset",x))
		}
	} else {
		if(!is.integer(x) || x > length(tsnames))
		{
			stop("Valid taxaset was not specified: ",x)
		} else {
			cmdobj$options = rbind(cmdobj$options,c("taxset",tsnames[x]))
		}
	}
	return(cmdobj)
}

addcmd.binary <- function(cmdobj,operand1,boolval)
{
	cmdobj$options = rbind(cmdobj$options,c(operand1,ifelse(boolval,"yes","no")  ))
	return(cmdobj)
}

addcmd.literal <- function(cmdobj,operand1,litval)
{
	cmdobj$options = rbind(cmdobj$options,c(operand1,litval))
	return(cmdobj)
}


# Non-censored rate test:
# arguments are options that can be specified:
#
# ... arguments for the commands(x,...) call.  add=? index=? replace=?
# 		allows the user to specify where this command should go.
#
addNonCensored <- function(obj,
							file=NULL,
							taxset="ALL",
							treeloop=FALSE,
							charloop=FALSE,
							fileappend=FALSE,
							filereplace=FALSE,
							usetempfile=FALSE,
							...)
{
	
	# convert to list for compatibility:
	if(!is.list(obj)){
		obj = list(obj)
	}
	
	if(hasCommands(obj[[1]]))
	{
		cmds = commands(obj[[1]])
		cat("New command being added after current ones:\n")
		for(jj in seq(length(cmds)))
			cat(cmds[jj],"\n")
		
	}
	
	# setup new command object:
	#
	newcmd = list(command=NULL,options=matrix(NA,ncol=2,nrow=0))
	newcmd$command = "cont"
	
	# add options
	newcmd = addcmd.taxaset(obj[[1]],newcmd,taxset) # taxaset
	newcmd = addcmd.binary(newcmd,"treeloop",treeloop) # treeloop
	newcmd = addcmd.binary(newcmd,"charloop",charloop) # charloop
	if(!usetempfile && !is.null(file) && is.character(file)) # file, append, replace
	{
		if(!file.exists(file))
			stop("Specified file,",file,", does not exist.")
		
		newcmd = addcmd.literal(newcmd,"file",file[1])
		newcmd = addcmd.binary(newcmd,"Append",fileappend)
		newcmd = addcmd.binary(newcmd,"Replace",filereplace)
		
	} else {
		if(usetempfile)
		{
			newcmd = addcmd.literal(newcmd,"file",tempfile())
			newcmd = addcmd.binary(newcmd,"Append",FALSE)
			newcmd = addcmd.binary(newcmd,"Replace",TRUE)		
		} else {
			newcmd = addcmd.literal(newcmd,"file","")
		}
	}

	# add commands to all members of the list:
	#
	for(ii in seq(length(obj)))
	{
		commands(obj[[ii]],...) <- write.brownie.string(newcmd)
	}
	
	if(length(obj) == 1)
		obj = obj[[1]]
	
	return(obj)
}



##### Censored Rate test:
#
addCensored <- function(obj,
						file=NULL,
						reps=0,
						taxset="ALL",
						treeloop=FALSE,
						charloop=FALSE,
						quiet=TRUE,
						usetempfile=FALSE,
						...)
{
	###### initialize
	if(!is.list(obj)){
		obj = list(obj)
	}							
	
	###### newcmd
	newcmd = list(command=NULL,options=matrix(NA,ncol=2,nrow=0))
	newcmd$command = "ratetest"
	
	newcmd = addcmd.taxaset(obj[[1]],newcmd,taxset) # taxaset
	newcmd = addcmd.literal(newcmd,"reps",as.character(reps))
	newcmd = addcmd.binary(newcmd,"treeloop",treeloop) # treeloop
	newcmd = addcmd.binary(newcmd,"charloop",charloop) # charloop
	newcmd = addcmd.binary(newcmd,"quiet",quiet) # charloop
	if(!usetempfile && !is.null(file) && is.character(file)) # file, append, replace
	{
		if(!file.exists(file))
			stop("Specified file,",file,", does not exist.")
		
		newcmd = addcmd.literal(newcmd,"file",file[1])		
	} else {
		if(usetempfile)
		{
			newcmd = addcmd.literal(newcmd,"file",tempfile())
		} else {
			newcmd = addcmd.literal(newcmd,"file","")
		}
	}
	
	
	# TODO: -Check for all trees whether taxasets are monophyletic or paraphyletic.
	#		 warn if some trees are not either.
	if(!all(sapply(obj,areTaxaMono,taxset)))
		warning("In some trees the taxaset,",taxset,", does not form a monophyletic clade.  Might want to remove those trees before running this analysis")
	
	
	# add commands to all members of the list:
	#
	for(ii in seq(length(obj)))
	{
		commands(obj[[ii]],...) <- write.brownie.string(newcmd)
	}
	
	###### return
	if(length(obj) == 1)
		obj = obj[[1]]
	
	return(obj)	
}




