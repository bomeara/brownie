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
	cmdobj
}

addcmd.binary <- function(cmdobj,operand1,boolval)
{
	cmdobj$options = rbind(cmdobj$options,c(operand1,ifelse(boolval,"yes","no")  ))
	cmdobj
}

addcmd.literal <- function(cmdobj,operand1,litval)
{
	cmdobj$options = rbind(cmdobj$options,c(operand1,litval))
	cmdobj
}

addcmd.model.discrete <- function(cmdobj,operand2,ratematrix,usingdata)
{
	# check if operand is valid and ratematrix makes sense
	if(!checkval.model.discrete(operand2))
		stop("Discete model specified is not a valid model: ",operand2)
	
	cmdobj$options = rbind(cmdobj$options,c("model",operand2))
	if(toupper(operand2) == "USER")
	{
		if(is.null(ratematrix))
			stop("Cannot set discrete model to USER without also specifying a rate matrix")
		
		if(!checkval.ratemat(ratematrix,usingdata))
			stop("Rate matrix: ",ratematrix," is not valid with data (levels=",levels(usingdata),").  See write.brownie.matrix for how to construct a rate matrix.")
		
		cmdobj$options = rbind(cmdobj$options,c("ratemat",ratematrix))
	}
	cmdobj
}

addcmd.model.continuous <- function(cmdobj,operand2)
{
	if(!checkval.model.continuous(operand2))
		stop("Continuos model specified is not a valid model: ",operand2)
	
	cmdobj$options = rbind(cmdobj$options,c("type",operand2))
	cmdobj
}

addcmd.freq <- function(cmdobj,operand2,statevect,usingdata)
{
	# check if operand is valid and state vector makes sense
	if(!check.freq(operand2))
		stop("Frequency model specified is not valid mode: ",operand2)
	
	cmdobj$options = rbind(cmdobj$options,c("freq",operand2))
	if(toupper(operand2) == "SET")
	{
		if(!check.statevector(statevect,usingdata))
			stop("Frequency vector,",statevect,", is not valid.")
		
		cmdobj$options = rbind(cmdobj$options,c("statevector",statevect))
	}
	cmdobj
}

addcmd.states <- function(cmdobj,statevect,usingdata)
{
	if(!check.statevector(statevect,usingdata))
		stop("State vector,",statevect,", is not valid.")
	
	cmdobj$options = rbind(cmdobj$options,c("states",statevect))
	
	cmdobj
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
		#if(!file.exists(file))
		#	stop("Specified file,",file,", does not exist.")
		
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
			newcmd = addcmd.literal(newcmd,"file","''")
		}
	}

	# add commands to all members of the list:
	#
	for(ii in seq(length(obj)))
	{
		commands(obj[[ii]],add=T,...) <- write.brownie.string(newcmd)
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
	
	if(length(taxset) != 0)
		for(jj in seq(length(taxset)))
			newcmd = addcmd.taxaset(obj[[1]],newcmd,taxset[jj]) # taxaset
		
	newcmd = addcmd.literal(newcmd,"reps",as.character(reps))
	newcmd = addcmd.binary(newcmd,"treeloop",treeloop) # treeloop
	newcmd = addcmd.binary(newcmd,"charloop",charloop) # charloop
	newcmd = addcmd.binary(newcmd,"quiet",quiet) # charloop
	if(!usetempfile && !is.null(file) && is.character(file)) # file, append, replace
	{
		#if(!file.exists(file))
		#	stop("Specified file,",file,", does not exist.")
		
		newcmd = addcmd.literal(newcmd,"file",file[1])		
	} else {
		if(usetempfile)
		{
			newcmd = addcmd.literal(newcmd,"file",tempfile())
		} else {
			newcmd = addcmd.literal(newcmd,"file","''")
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
		commands(obj[[ii]],add=T,...) <- write.brownie.string(newcmd)
	}
	
	###### return
	if(length(obj) == 1)
		obj = obj[[1]]
	
	return(obj)	
}


##### Discrete character evolution
#
# NOTE: will this automatically use the first data column?
#
addDiscrete <- function(obj,
						file=NULL,
						model=brownie.models.discrete()[1],
						model.state=NULL,
						ratemat=NULL,
						freq=brownie.freqs()[1],
						statevector=NULL,
						treeloop=FALSE,
						charloop=FALSE,
						allchar=FALSE,
						variable=FALSE,
						reconstruct=FALSE,
						breaknum=0,
						fileappend=FALSE,
						filereplace=FALSE,
						globalstates=FALSE,
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
	
	# RATE MATRIX
	dat = NULL
	if(is.null(model.state)){
		colind = head(which(datatypes(obj[[1]]) == discData()),1)
		if(length(colind)==0)
			stop("Could not find data column to use.")
		dat = tdata(obj[[1]],'tip')[,colind,drop=T]
	} else {
		dat = tdata(obj[[1]],'tip')[,colind,drop=T]
	}
	newcmd = addcmd.discrete(newcmd,model,ratemat,dat)
	
	# FREQUENCIES:
	newcmd = addcmd.freq(newcmd,freq,statevector,dat)
	
	newcmd = addcmd.binary(newcmd,"treeloop",treeloop) # treeloop
	newcmd = addcmd.binary(newcmd,"charloop",charloop) # charloop
	newcmd = addcmd.binary(newcmd,"allchar",allchar) # allchar
	newcmd = addcmd.binary(newcmd,"variable",variable) # variable
	newcmd = addcmd.binary(newcmd,"reconstruct",reconstruct) # reconstruct
	newcmd = addcmd.binary(newcmd,"breaknum",as.character(breaknum)) # break number
	
	# FILE:
	if(!usetempfile && !is.null(file) && is.character(file)) # file, append, replace
	{
		#if(!file.exists(file))
		#	stop("Specified file,",file,", does not exist.")
		
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
	newcmd = addcmd.binary(newcmd,"GlobalStates",globalstates) # allchar
		
	# add commands to all members of the list:
	#
	for(ii in seq(length(obj)))
	{
		commands(obj[[ii]],add=T,...) <- write.brownie.string(newcmd)
	}
	
	###### return
	if(length(obj) == 1)
		obj = obj[[1]]

	return(obj)
}


##### Start log file recording
#
addStartLog <- function(obj,
							file=NULL,
						fileappend=FALSE,
						filereplace=FALSE,
						usetempfile=FALSE,
						...)
{
	###### initialize
	if(!is.list(obj)){
		obj = list(obj)
	}
	
	newcmd = list(command=NULL,options=matrix(NA,ncol=2,nrow=0))
	newcmd$command = "log start"
	
	# FILE:
	if(!usetempfile && !is.null(file) && is.character(file)) # file, append, replace
	{
		#if(!file.exists(file))
		#	stop("Specified file,",file,", does not exist.")
		
		newcmd = addcmd.literal(newcmd,"file",file[1])
		newcmd = addcmd.binary(newcmd,"Append",fileappend)
		newcmd = addcmd.binary(newcmd,"Replace",filereplace)
		
	} else {
		if(usetempfile)
		{
			newcmd = addcmd.literal(newcmd,"file",tempfile())
			newcmd = addcmd.binary(newcmd,"Append",FALSE)
			newcmd = addcmd.binary(newcmd,"Replace",TRUE)		
		}
	}	
	
	for(ii in seq(length(obj)))
	{
		if(nrow(newcmd$options)!=0){
			commands(obj[[ii]],add=T,...) <- write.brownie.string(newcmd)
		} else {
			warning("Could not find file,",file," to log to.")
		}
	}
		
	###### return
	if(length(obj) == 1)
		obj = obj[[1]]

	return(obj)
}

##### End log file recording
#
addEndLog <- function(obj)
{
	###### initialize
	if(!is.list(obj)){
		obj = list(obj)
	}
	

	for(ii in seq(length(obj)))
	{
			commands(obj[[ii]],add=T,...) <- "log stop;"
	}
		
	###### return
	if(length(obj) == 1)
		obj = obj[[1]]

	return(obj)	
}


##### Continuous models
#
addModel <- function(obj,
					type=brownie.models.continuous()[1],
					states=NULL, # default to all different
					changes=NULL, # defaults to 10
					model.state=NULL)
{
	###### initialize
	if(!is.list(obj)){
		obj = list(obj)
	}
	
	
	newcmd = list(command=NULL,options=matrix(NA,ncol=2,nrow=0))
	newcmd$command = "model"
	
	newcmd = addcmd.model.continuous(newcmd,type)
	
	if(!is.null(states))
	{
		dat = NULL
		if(is.null(model.state)){
			colind = head(which(datatypes(obj[[1]]) == discData()),1)
			if(length(colind)==0)
				stop("Could not find data column to use.")
			dat = tdata(obj[[1]],'tip')[,colind,drop=T]
		} else {
			dat = tdata(obj[[1]],'tip')[,colind,drop=T]
		}
		
		newcmd = addcmd.states(newcmd,states,dat)
	}
	
	if(!is.null(changes) && !is.na(as.integer(changes)) && as.integer(changes) > 0)
		newcmd = addcmd.literal(newcmd,"changes",as.character(changes))
	
	
	for(ii in seq(length(obj)))
	{
		commands(obj[[ii]],add=T) <- write.brownie.string(newcmd)
	}
		
	###### return
	if(length(obj) == 1)
		obj = obj[[1]]

	return(obj)	
}


