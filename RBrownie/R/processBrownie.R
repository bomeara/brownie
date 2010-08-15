#--------------------------------------------------
# Supplemental brownie functions:
# -	Need to load these things first
#
#--------------------------------------------------

## Valid 'Datatypes':
# (for use with @datatypes slot)

# datatypes:
discData <- function() { return("discrete") }
contData <- function() { return("cont") }
taxaData <- function() { return("taxset") }
genericData <- function() { return("undef") }

brownie.datatypes <-function()
{
	brownie_datatype_options = c(taxaData(),contData(),discData(),genericData())
	return(brownie_datatype_options)
}

# convert between datatypes:
# helper functions for addData
as.discData <- function(dat)
{
	if(!is.factor(dat))
		dat = as.factor(dat)
	return(dat)
}

is.discData <- function(dat)
{
	return(is.factor(dat))
}

as.contData <- function(dat)
{
	if(!is.numeric(dat))
	{
		if(is.factor(dat))
		{
			dat = as.numeric(levels(dat)[as.numeric(dat)])
		} else {
			# assume there are characters:
			dat = as.numeric(dat)
		}
	}
	return(dat)
}

is.contData <- function(dat)
{
	return(is.numeric(dat))
}




## Valid Brownie Commands:
#
#



## Valid brownie models
#
brownie.models.continuous <- function(with.desc=F)
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


brownie.models.discrete <- function(with.desc=F)
{
	brownie_model_options = c(
		"EQUAL",
		"REV",
		"NONREV",
		"USER"
	)
	
	brownie_model_description = c(
		"Model where all rates are equal",
		"Reversible model where q_ij == q_ji for all i and j",
		"Non-reversible model where all rates, q_ij, can vary independently.",
		"User-specified model."
	)
	
	if(with.desc)
		brownie_model_options = cbind(brownie_model_options,brownie_model_description)
	
	return(brownie_model_options)
}



brownie.freqs <- function(with.desc=F)
{
	brownie_freqs_options = c(
		"EMPIRICAL",
		"EQUILIBRIUM",
		"UNIFORM",
		"OPTIMIZE",
		"SET"
	)
	
	brownie_freqs_description = c(
		"Frequencies based on empirical distribution of states at the tips",	
		"The frequencies expected with the optimized rate matrix given infinitely long branches",
		"Frequencies all the same",
		"Frequencies are optimized as part of the model",
		"User-specified frequencies"
	)
	
	if(with.desc)
		brownie_freqs_options = cbind(brownie_freqs_options,brownie_freqs_description)
	
	return(brownie_freqs_options)
}

# tipvariance types
brownie.tvtypes <- function(with.desc=F)
{
	brownie_freqs_options = c(
		"None",
		"Given",
		"Same"
	)
	
	brownie_freqs_description = c(
		"Assume no tipvariance",	
		"User-given tip variances for each taxon",
		"Estimate one tipvariance across all taxa."
	)
	
	if(with.desc)
		brownie_freqs_options = cbind(brownie_freqs_options,brownie_freqs_description)
	
	return(brownie_freqs_options)
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

checkval.model.discrete <- function(optstr)
{
	# TODO: figure out if brownie is in fact case-insensitive like I'm assuming
	#
	return( length(grep(tolower(optstr),tolower(brownie.models.discrete())))>0 )
}

checkval.model.continuous <- function(optstr)
{
	# TODO: figure out if brownie is in fact case-insensitive like I'm assuming
	#
	return( length(grep(tolower(optstr), tolower(brownie.models.continuous())))>0 )
}

checkval.ratemat <- function(optstr,factvect)
{
	# rate matrix needs to be in format: (a b c d) (no diagonals)
	levs = levels(factvect)
	numberneeded = (length(levs)*length(levs)) - length(levs)
	if(numberneeded > 2)
		return(FALSE)
	
	inside = sub("^\\((.*)\\)$","\\1",optstr)
	tokens = strsplit(inside," ")[[1]]
	
	# match these:
	validtokens = length( grep("^([0-9A-Za-z]|[-+]?[0-9]\\d{0,2}(\\.\\d{1,10})?%?)$",tokens) )
	
	return( (validtokens == numberneeded) )
	
}

check.freq <- function(optstr)
{
	return( length(grep(tolower(optstr),tolower(brownie.freqs())))>0 )
}

check.statevector <- function(optstr,factvect)
{
	levs = levels(factvect)
	inside = sub("^\\((.*)\\)$","\\1",optstr)
	tokens = strsplit(inside," ")[[1]]
	validtokens = tokens[grep("^([-+]?[0-9]\\d{0,2}(\\.\\d{1,10})?%?)$",tokens)]
	
	return( (length(levs) == length(validtokens)) )
}



#
#------------------------------------------------------------------
# Parsing Brownie blocks into tokens for validation.
#------------------------------------------------------------------
#
# does the string contain char
haschar<-function(str,char)
{
	as.logical(length(grep(char,str)))
}

# roughly, is the char in the back or front
isback <- function(str,char)
{
	as.logical( round (which(strsplit(str,"")[[1]] == char ) / nchar(str)) )
}


# assume string is "A<char>B" format, with A or B optional, but not both
get.left.right <- function(str,char)
{
	outvals = strsplit(str,char)[[1]]
	if(length(outvals) > 2)
		stop("Problem parsing command: more than one operator in a token")
	
	if(length(outvals)==1)
	{
		if(isback(str,char))
			return( c(outvals,"") )
		
		return( c("",outvals) )
	}
	
	return( outvals )
}

check.empty <- function(bb) sapply(bb,function(i) i=="")

rm.ending <- function(str)
{
	if(haschar(str,";$"))
		 str = sub(";$","",str)
	
	return(str)
}


read.brownie.string <- function(txtstr,operator="=")
{
	command = ""
	opts = matrix(character(0),ncol=2,nrow=0)
	optcount = 0

	tokens = strsplit(txtstr,"\\s")[[1]]
	command = head(tokens,1)
	tokens = tokens[-1]

	operator="=" # should be const
	isleft=TRUE # look at left of the equals?
	for(token in tokens)
	{
		#print(token)
		if(token == operator)
		{
			if(isleft)
				stop("Problem parsing command: '=' before first operand")
			next
		}
	
		if(isleft)
		{
			if(haschar(token,operator))
			{
				# process first operand
				both = get.left.right(token,operator)
				both.isempty = check.empty(both)
				if(all(both.isempty))
					warning("Empty token!")
				
				if(both.isempty[1])
					stop("Problem parsing command: 1st operand is empty")
				
				opts = rbind(opts,both)
				if(both.isempty[2])
					isleft = FALSE
				
			} else {
				opts = rbind(opts,c(token,""))
				isleft = FALSE
			}
			optcount = optcount + 1
			
		} else {
			
			if(haschar(token,operator))
			{
				# process first operand
				both = get.left.right(token,operator)
				both.isempty = check.empty(both)
				if(all(both.isempty))
					warning("Empty token!")

				if(both.isempty[2])
					stop("Problem parsing command: 2nd operand is empty")
				
				opts[optcount,2] = both[2]
			} else {
				opts[optcount,2] = token
			}
		
			isleft = TRUE
		}
	}
	opts = unname(opts)
	opts[,2] = unname(sapply(opts[,2],rm.ending))
	colnames(opts) <- c("operand.1","operand.2")
	
	return(list(command=command,options=opts))
}


# assume obj is list with $commands and $options
write.brownie.string <- function(obj,operator="=")
{
	if(is.null(obj$command))
		stop("Couldn't find 'command' slot for input")
	
	retstr = obj$command
	
	if(!is.null(obj$options))
	{
		tmp=apply(obj$options,1,paste,collapse="=")
		tmp = paste(tmp,collapse=" ")
		retstr = paste(retstr,tmp)
	}

	return( paste(retstr,";",sep="") )
}



#--------------------------------------------------------------
# Helper functions to write brownie cmd strings.
#--------------------------------------------------------------

# mat should be in format:
# to -> 0	1	2
# from	---------
# 	| | -	a	0.5
#   v | b	-	c
# 	  | 0.0	a	-
#
# Example:
# --------
# mata = matrix(c("-","a","0.5","b","-","c","0.0","a","-"),nrow=3,byrow=T)
# write.brownie.matrix(mata)
#
write.brownie.matrix<-function(mat)
{
	retstr = character(0)
	for(ii in seq(nrow(mat)))
		for(jj in seq(ncol(mat)))
			if(ii!=jj)
				retstr = append(retstr,mat[ii,jj])
	
	retstr = paste(retstr,collapse=" ")
	retstr = sprintf("(%s)",retstr)
	return(retstr)
}

# brownie vector:
# Example:
# -------
# write.brownie.vector(c(0.4,0.6))
#
write.brownie.vector <- function(vect,rescale=T)
{
	if(rescale)
		vect = vect / sum(vect)
	
	retstr = sprintf("(%s)",paste(vect,collapse=" "))
	return(retstr)
}


