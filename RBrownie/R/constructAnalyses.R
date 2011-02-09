#
#--------------------------------------------------
# Construct Phylogenetic analyses for RBrownie
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
run.analysis <- function(filename)
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
	
	return( .Call("doAnalysis",filename,PACKAGE="RBrownie") )
}

# Run a brownie file as it is.  If you just want to see if brownie can
# execute the file (without crashing), set debugging to TRUE.
#
run.asis <- function(brobj,brfile=NULL,logfile=NULL,debugging=F)
{
	if(!debugging && !hasCommands(brobj))
	{
		warning("This brownie object has no commands")
		return(list(textout=character(0),treesout=character(0)))
	}
	if(is.null(brfile))
		brfile=tempfile()
	
	if(!is.null(logfile))
	{
		brobj = addStartLog(brobj,file=logfile,filereplace=T)
		brobj = addEndLog(brobj)
	}
	
	writeBrownie(brobj,brfile)
	
	outtext = run.analysis(brfile)
	
	return(outtext)
}



# Do ancestral state reconstruction
#
# @param return phylo4d_ext format or non-extended phylo4d
# @file file to dump output to
# @... options to be passed to addDiscrete...
#
runDiscrete <- function(brobj,outfile=NULL,brfile=NULL,models=brownie.models.discrete()[1],freqs=brownie.freqs()[1],...)
{
	
	# sanity checks:
	if(length(models) != length(freqs))
		stop("Must specify the same number of models and frequencies\nType: discreteOptions() for more information.")
			
	# lists are easier to work with
	if(!is(brobj,'list'))
	{
		brobj = list(brobj)
	}
	# convert to brownie if needed
	if(!is(brobj[[1]],"brownie"))
	{
		brobj = brownie(brobj)
	}
	# clear out all other commands:
	brobj = clearCommands(brobj)
	
	cmds = character(seq(length(models)))
	for(ii in seq(length(models)))
	{
		brobj = addDiscrete(brobj,file=outfile,model=models[ii],freq=freqs[ii],...)
		cmds[ii] = tail(commands(brobj),1)
	}
	
	if(is.null(brfile))
		brfile=tempfile()
	
	writeBrownie(brobj,brfile)
	
	outtext = run.analysis(brfile)
	
	
	rettxt = read.discrete.output(txt=scan.textout(outtext$textout))
	rettr = scan.treesout(outtext$treesout)
	
	retlist = list(stats=rettxt,trees=rettr,cmds=cmds)
	
	return( retlist )
}


# Do uncensored rate test ("continuous","opt")
#
runNonCensored <- function(brobj,outfile=NULL,brfile=NULL,
							models=brownie.models.continuous()[1],
							taxsets=character(0),
							states=character(0),
							changes=character(0),...)
{
	
	# sanity checks:
	if(length(taxsets)>1 && length(models)!=length(taxsets))
		stop("If you do specify taxasets, the vector must be the same length as the models")
	
	cchanges = NULL
	cstates = NULL
	if(length(changes)==1)
		cchanges = rep(changes,length(models))
	
	if(length(states)==1)
		cstates = rep(states, length(models))
	
	# lists are easier to work with
	if(!is(brobj,'list'))
	{
		brobj = list(brobj)
	}
	# convert to brownie if needed
	if(!is(brobj[[1]],"brownie"))
	{
		brobj = brownie(brobj)
	}
	# clear out all other commands:
	brobj = clearCommands(brobj)
		
	for(ii in seq(length(models)))
	{
		copt = NULL
		sopt = NULL
		if(!is.null(cstates))
			sopt = cstates[ii]
		
		if(!is.null(cchanges))
			copt = cchanges[ii]
		
		brobj = addModel(brobj,type=models[ii],changes=copt,states=sopt)
		
		if(length(taxsets)==0){
			brobj = addNonCensored(brobj,file=outfile,...)			
		} else {
			brobj = addNonCensored(brobj,file=outfile,taxset=taxsets[ii],...)
		}
	}

	if(is.null(brfile))
		brfile=tempfile()
	
	writeBrownie(brobj,brfile)
	
	outtext = run.analysis(brfile)
	if(length(outtext$textout)>0)
	{
		# NOTE: this will be a problem if outdat == NULL (if the wrong model is used).
		outdat = read.continuous.output(txt=scan.textout(outtext$textout[1]))
		
		if(length(outtext$textout)>1){
			for(xx in seq(from=2,to=length(outtext$textout)))
			{
				tmpdf = read.continuous.output(txt=scan.textout(outtext$textout[xx]))
				if( !(ncol(tmpdf)==0 || nrow(tmpdf)==0) )
					outdat = merge(outdat,tmpdf,all=T)
			}
		}
	}
	
	if(!is.null(outdat$Model))
		outdat = outdat[order(outdat$Model),]
	
	return( outdat )	
}


# Do censored rate test ("ratetest")
# @param ellipses are parameters to be passed to addCensored
# 
#
runCensored <- function(brobj,brfile=NULL,file=NULL,...)
{	
	
	# sanity checks
	if(!is.null(file))
		if(file.exists(file))
			stop("File ",file," already exists and will cause this test to fail to run")
	
	# lists are easier to work with
	if(!is(brobj,'list'))
	{
		brobj = list(brobj)
	}
	# convert to brownie if needed
	if(!is(brobj[[1]],"brownie"))
	{
		brobj = brownie(brobj)
	}
	# clear out all other commands:
	brobj = clearCommands(brobj)
	
	# only need one
	brobj = addCensored(brobj,file=file,...)

	if(is.null(brfile))
		brfile=tempfile()
	
	writeBrownie(brobj,brfile)	
	outtext = run.analysis(brfile)
	
	outdat = read.ratetest.output(txt=scan.textout(outtext$textout))
	
	return( outdat )
	
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


addcmd.tvtype <- function(cmdobj,type)
{
	if(!tolower(type) %in% tolower(brownie.tvtypes()))
		stop("tipvariance type,",type,", is not a valid tipvariance.")
	
	cmdobj$options = rbind(cmdobj$options,c("type",type))
	
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
							taxset=character(0),
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
		#cat("New command being added after current ones:\n")
		#for(jj in seq(length(cmds)))
		#	cat(cmds[jj],"\n")
	}
	
	# setup new command object:
	#
	newcmd = list(command=NULL,options=matrix(NA,ncol=2,nrow=0))
	newcmd$command = "cont"
	
	# add options
	if(length(taxset)!=0 && is.character(taxset))
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
			#newcmd = addcmd.literal(newcmd,"file","''")
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
						taxset=character(0),
						treeloop=FALSE,
						charloop=FALSE,
						quiet=FALSE,
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
	
	if(length(taxset) != 0 && is.character(taxset))
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
			#newcmd = addcmd.literal(newcmd,"file","''")
		}
	}
	
	
	# TODO: -Check for all trees whether taxasets are monophyletic or paraphyletic.
	#		 warn if some trees are not either.
	if(length(taxset) != 0)
		for(jj in seq(length(taxset)))
			if(!all(sapply(obj,areTaxaMono,taxset[jj])))
				warning("In some trees the taxaset,",taxset[jj],", does not form a monophyletic clade.  Might want to remove those trees before running this analysis")
	
	
	# add commands to all members of the list:
	#
	for(ii in seq(length(obj)))
	{
		commands(obj[[ii]],add=T) <- write.brownie.string(newcmd)
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
	newcmd$command = "discrete"
	
	# RATE MATRIX
	dat = NULL
	if(is.null(model.state)){
		colind = head(which(datatypes(obj[[1]]) == discData()),1)
		if(length(colind)==0)
			stop("Could not find data column to use.")
		dat = tdata(obj[[1]],'tip')[,colind,drop=T]
	} else {
		if(is.integer(model.state)){
			colind=model.state
		} else {
			colind = which(colnames(tdata(obj[[1]],'tip')) == model.state)
		}
		dat = tdata(obj[[1]],'tip')[,colind,drop=T]
	}
	newcmd = addcmd.model.discrete(newcmd,model,ratemat,dat)
	
	# FREQUENCIES:
	newcmd = addcmd.freq(newcmd,freq,statevector,dat)
	
	newcmd = addcmd.binary(newcmd,"treeloop",treeloop) # treeloop
	newcmd = addcmd.binary(newcmd,"charloop",charloop) # charloop
	newcmd = addcmd.binary(newcmd,"allchar",allchar) # allchar
	newcmd = addcmd.binary(newcmd,"variable",variable) # variable
	newcmd = addcmd.binary(newcmd,"reconstruct",reconstruct) # reconstruct
	newcmd = addcmd.literal(newcmd,"breaknum",as.character(breaknum)) # break number
	
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
			#newcmd = addcmd.literal(newcmd,"file","")
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
	newcmd$command = "log"
	
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
	return( addLiteral(obj,"log stop;") )
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
			if(is.integer(model.state)){
				colind=model.state
			} else {
				colind = which(colnames(tdata(obj[[1]],'tip')) == model.state)
			}
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


# add opt command:
#	" Returns the likelihood and AICc under the current model"
#
addOpt <- function(obj,
					file=NULL,
					taxset=character(0),
					treeloop=FALSE,
					charloop=FALSE,
					usetempfile=FALSE)
{
	###### initialize
	if(!is.list(obj)){
		obj = list(obj)
	}							
	
	###### newcmd
	newcmd = list(command=NULL,options=matrix(NA,ncol=2,nrow=0))
	newcmd$command = "opt"
	
	if(length(taxset) != 0 && is.character(taxset))
		for(jj in seq(length(taxset)))
			newcmd = addcmd.taxaset(obj[[1]],newcmd,taxset[jj]) # taxaset
		
	newcmd = addcmd.binary(newcmd,"treeloop",treeloop) # treeloop
	newcmd = addcmd.binary(newcmd,"charloop",charloop) # charloop
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
			#newcmd = addcmd.literal(newcmd,"file","''")
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
		commands(obj[[ii]],add=T) <- write.brownie.string(newcmd)
	}
	
	###### return
	if(length(obj) == 1)
		obj = obj[[1]]
	
	return(obj)	
}


addTipvariance <- function(obj,type=brownie.tvtypes()[1])
{
	###### initialize
	if(!is.list(obj)){
		obj = list(obj)
	}							
	
	###### newcmd
	newcmd = list(command=NULL,options=matrix(NA,ncol=2,nrow=0))
	newcmd$command = "tipvariance"
	
	newcmd = addcmd.tvtype(newcmd,type) # treeloop
	
	# add commands to all members of the list:
	#
	for(ii in seq(length(obj)))
	{
		commands(obj[[ii]],add=T) <- write.brownie.string(newcmd)
	}
	
	###### return
	if(length(obj) == 1)
		obj = obj[[1]]
	
	return(obj)		
}


# add 'choose' command
# from brownie: message="Usage: Choose [tree=<integer>] [char=<integer>]\n\n";
#
addChoose <- function(obj,
					char=NULL,
					tree=NULL,
					discchar=NULL)
{
	###### initialize
	if(!is.list(obj)){
		obj = list(obj)
	}
	
	###### newcmd
	newcmd = list(command=NULL,options=matrix(NA,ncol=2,nrow=0))
	newcmd$command = "choose"
	
	if(!is.null(char))
		newcmd = addcmd.literal(newcmd,"char",as.character(char))
	
	if(!is.null(tree))
		newcmd = addcmd.literal(newcmd,"tree",as.character(tree))
	
	if(!is.null(discchar))
		newcmd = addcmd.literal(newcmd,"d",as.character(discchar))
	
	# add commands to all members of the list:
	#
	for(ii in seq(length(obj)))
	{
		commands(obj[[ii]],add=T) <- write.brownie.string(newcmd)
	}
	
	###### return
	if(length(obj) == 1)
		obj = obj[[1]]
	
	return(obj)		
}


addNumopt <- function(obj,
					iter=NULL,
					toler=NULL,
					randstart=NULL,
					seed=NULL,
					stepsize=NULL,
					detail=FALSE,
					redo=FALSE,
					giveupfactor=NULL
)
{
	###### initialize
	if(!is.list(obj)){
		obj = list(obj)
	}
	
	###### newcmd
	newcmd = list(command=NULL,options=matrix(NA,ncol=2,nrow=0))
	newcmd$command = "numopt"
	
	if(!is.null(iter) && is.numeric(iter))
		newcmd = addcmd.literal(newcmd,"iter",as.character(iter))
	
	if(!is.null(toler) && is.numeric(toler))
		newcmd = addcmd.literal(newcmd,"toler",as.character(toler))

	if(!is.null(randstart) && is.numeric(randstart))
		newcmd = addcmd.literal(newcmd,"randstart",as.character(randstart))
	
	if(!is.null(seed) && is.numeric(seed))
		newcmd = addcmd.literal(newcmd,"seed",as.character(seed))
	
	if(!is.null(stepsize) && is.numeric(stepsize))
		newcmd = addcmd.literal(newcmd,"stepsize",as.character(stepsize))
	
	newcmd = addcmd.binary(newcmd,"detail",detail)	
	newcmd = addcmd.binary(newcmd,"redo",redo)

	if(!is.null(giveupfactor) && is.numeric(giveupfactor))
		newcmd = addcmd.literal(newcmd,"giveupfactor",as.character(giveupfactor))
	
	
	# add commands to all members of the list:
	#
	for(ii in seq(length(obj)))
	{
		commands(obj[[ii]],add=T) <- write.brownie.string(newcmd)
	}
	
	###### return
	if(length(obj) == 1)
		obj = obj[[1]]
	
	return(obj)	
}

# Add heuristic search.
#
addHS <- function(obj,
					nreps=NULL,#=5,
					rearrlimit=NULL,
					maxnumspp=NULL,#=100,
					minnumspp=NULL,#=1,
					minsamp=NULL,#=3,
					structwt=NULL,#=0.5,
					pthreshold=NULL,#=1.0,
					subsample=NULL,#=2.0,
					movefreq=NULL,#="(0.8 0.01 0.01 0.01 0.17 0.0)",
					file=NULL,#="besttrees.tre",
					showtries=NULL,#=No,
					coal=NULL,#=No,
					ms=NULL,#=No,
					aic.mode=NULL,#=4,
					gridwidth=NULL,#=10.0,
					gridsize=NULL,#=15.0,
					maxrecursions=NULL,#=20,
					branch.export=NULL,#=2
					usetempfile=FALSE
)
{
	###### initialize
	if(!is.list(obj)){
		obj = list(obj)
	}
	
	###### newcmd
	newcmd = list(command=NULL,options=matrix(NA,ncol=2,nrow=0))
	newcmd$command = "hs"
	
	if(!is.null(nreps) && is.numeric(nreps))
		newcmd = addcmd.literal(newcmd,"NReps",as.character(nreps))
	
	if(!is.null(rearrlimit) && is.numeric(rearrlimit))
		newcmd = addcmd.literal(newcmd,"RearrLimit",as.character(rearrlimit))

	if(!is.null(maxnumspp) && is.numeric(maxnumspp))
		newcmd = addcmd.literal(newcmd,"MaxNumSpp",as.character(maxnumspp))
	
	if(!is.null(minnumspp) && is.numeric(minnumspp))
		newcmd = addcmd.literal(newcmd,"MinNumSpp",as.character(minnumspp))
		
	if(!is.null(minsamp) && is.numeric(minsamp))
		newcmd = addcmd.literal(newcmd,"MinSamp",as.character(minsamp))
	
	if(!is.null(structwt) && is.numeric(structwt))
		newcmd = addcmd.literal(newcmd,"StructWt",as.character(structwt))
		
	if(!is.null(pthreshold) && is.numeric(pthreshold))
		newcmd = addcmd.literal(newcmd,"PThreshold",as.character(pthreshold))
	
	if(!is.null(subsample) && is.numeric(subsample))
		newcmd = addcmd.literal(newcmd,"SubSample",as.character(subsample))

	if(!is.null(movefreq) && (check.statevector(movefreq,factor(1:5)) || check.statevector(movefreq,factor(1:6))) )
		newcmd = addcmd.literal(newcmd,"MoveFreq",as.character(movefreq))
		
	if(!usetempfile && !is.null(file) && is.character(file)) # file, append, replace
	{
		newcmd = addcmd.literal(newcmd,"file",file[1])
	} else {
		if(usetempfile)
			newcmd = addcmd.literal(newcmd,"file",tempfile())
	}
	
	if(!is.null(showtries) && is.logical(showtries))
		newcmd = addcmd.binary(newcmd,"showtries",showtries)

	# Be careful when using this option: it relies on external program
	if(!is.null(ms) && is.logical(ms))
		newcmd = addcmd.binary(newcmd,"MS:",ms)

	if(!is.null(coal) && is.logical(coal))
		newcmd = addcmd.binary(newcmd,"COAL:",coal)
										
	if(!is.null(aic.mode) && (aic.mode %in% seq(0,4)))
		newcmd = addcmd.literal(newcmd,"AIC_mode",as.character(aic.mode))

	if(!is.null(gridwidth) && is.numeric(gridwidth))
		newcmd = addcmd.literal(newcmd,"GridWidth",as.character(gridwidth))
		
	if(!is.null(gridsize) && is.numeric(gridsize))
		newcmd = addcmd.literal(newcmd,"GridSize",as.character(gridsize))
		
	if(!is.null(maxrecursions) && is.numeric(maxrecursions))
		newcmd = addcmd.literal(newcmd,"MaxRecursions",as.character(maxrecursions))
		
	if(!is.null(branch.export) && is.numeric(branch.export))
		newcmd = addcmd.literal(newcmd,"Branch_export",as.character(branch.export))
	
	
	for(ii in seq(length(obj)))
	{
		commands(obj[[ii]],add=T) <- write.brownie.string(newcmd)
	}
	
	###### return
	if(length(obj) == 1)
		obj = obj[[1]]
	
	return(obj)	
}


# This adds the cmdstr to the commands slot.  It does no checking.
# 
addLiteral <- function(obj,cmdstr=NULL)
{
	###### initialize
	if(!is.list(obj)){
		obj = list(obj)
	}
	
	# add ending semicolon
	if( length(grep(";$",cmdstr))==0 )
	{
		cmdstr=paste(cmdstr,";",sep="")
	}
	
	# add commands to all members of the list:
	#
	for(ii in seq(length(obj)))
	{
		commands(obj[[ii]],add=T) <- cmdstr
	}
	
	###### return
	if(length(obj) == 1)
		obj = obj[[1]]
	
	return(obj)
}


