#---------------------------------------------
# Process Nexus files  
# -	Extra methods for extracting different information
# 	from nexus-formatted files.
#
#---------------------------------------------


# Strip comments from a tree string
#
.strip.tree.comments <- function(text=NULL)
{
	
}


# Method to read the first comment in a line in the format '[&...]'
# This is for reading tree weights chiefly
#
# Example:
# get.nexus.comments("example.txt")->lala
#
get.nexus.comments<-function(finput,text=NULL)
{
	
	if(!is.null(text)){
		# TODO: check text for newlines and split on them if they exist.
		rawtext=text
	} else {
		if(!file.exists(finput))
			stop("Assuming finput is a file and could not find it")
		
		rawtext = scan(finput,what=character(0),strip.white=T,sep="\n")		
	}
	
	# TODO: return named pair {treename, comment}
	comments = character(0)
	for(ii in seq(length(rawtext)))
	{
		junk =  gsub("^.*?\\[(.*?)\\].*$","\\1",rawtext[ii])
		if(length( grep("^&(.*)$",junk) ) != 0)
			comments = append(comments,junk)
	}
	
	return(comments)	
}


has.weights <- function(finput,text=NULL)
{
	if(!is.null(text)){
		# TODO: check text for newlines and split on them if they exist.
		rawtext=text
	} else {
		if(!file.exists(finput))
			stop("Assuming finput is a file and could not find it")
		
		rawtext = scan(finput,what=character(0),strip.white=T,sep="\n")		
		rawtext = read.nexus.block(txt=rawtext,block="trees",silent=T)
	}
	
	comments = grep("(&lnP|&W)( |=)",rawtext)
	return( length(comments)>0 )
}


# get tree weights from file or tree string
# Assuming this format: 
# [.... &W -122235 ........]
#
get.tree.weights <- function(finput,text=NULL,starters=c("&W","&lnP"),splitchar="( |=)")
{
	if(!is.null(text)){
		# TODO: check text for newlines and split on them if they exist.
		rawtext=text
	} else {
		if(!file.exists(finput))
			stop("Assuming finput is a file and could not find it")
		
		rawtext = scan(finput,what=character(0),strip.white=T,sep="\n")		
		rawtext = read.nexus.block(txt=rawtext,block="trees",silent=T)
	}
	
	comments = get.nexus.comments(text=rawtext)
	tokens = strsplit(comments,splitchar)
	find.weight <- function(ii,ww)
	{
		return(ii[(which(toupper(ii) %in% toupper(ww)))+1])
	}
	weighs = unlist(lapply(tokens,find.weight,starters))
	
	return(as.numeric(weighs))
}



# Internal function - split tokens by a certain character
.split.tokens <- function(txt,char)
{
	nlines = length(txt)
	newtokens = character(nlines)
	charlines = unname(sapply(txt,function(i) length(grep(char,i))))
	curline=1
	for(ll in seq(nlines))
	{
		if(charlines[ll] == 0){
			newtokens[curline] = txt[ll]
			curline = curline + 1
		}else{
			tmp = strsplit(txt[ll],char)[[1]]
			for(tmpline in tmp){
				#if(tmpline!=""){
					#newtokens[curline] = paste(tmpline,char,sep="")
					newtokens[curline] = tmpline
					curline = curline + 1
				#}
			}
		}
	}
	return (newtokens)
}


# Internal function - get all content within a nexus block
.get.nexus.block.inds <- function(filename,blockname,text=NULL)
{
	# choose character vector
	if(!is.null(text))
	{
		filetext = text
	} else {
		filetext = scan(filename,what="character",sep="\n",strip.white=T)
	}
	filetext = tolower(filetext)
	
	search.for=paste(paste("begin",tolower(blockname),sep="\\s"),"[\\s]{0,};",sep="")
	start.ind = grep( search.for ,filetext,ignore.case=T)
	if(length(start.ind) == 0)
		return (integer(0))
	
	end.ind = grep("end;",filetext,ignore.case=T)
	end.ind = end.ind[head(which(end.ind > start.ind),1)]

	return( c((start.ind),(end.ind)) )
}



# Internal:
# Get text of a read nexus block
read.nexus.block<-function(finput,txt=NULL,block,rm.comments=F,silent=F)
{
	if(!is.null(txt)){
		# Using the text argument is not recommended
		rawtext=txt
	} else {
		
		if(!file.exists(finput))
			stop("Assuming finput is a file and could not find it")
		
		rawtext = scan(finput,what=character(0),strip.white=T,sep="\n",quiet=silent)		
	}
	
	inds = .get.nexus.block.inds(blockname=block,text=rawtext)
	if(length(inds)==0)
	{
		if(!silent) warning(paste("This file has",block, "no block"))
		return (character(0))
	}
	
	rawtext = rawtext[(inds[1]+1):(inds[2]-1)]
	# TODO: split up newlines if they exist
	
	if(rm.comments)
	{
		# assume comments are start a line with [ 
		# and end a line with ]
		#
		comment.starts = grep("^\\[",rawtext)
		comment.ending = grep("\\]$",rawtext)
		if(length(comment.starts) == length(comment.ending)){
			comment.pairs = cbind(comment.starts,comment.ending)
		} else {
			comment.pairs = cbind(character(0),character(0))
		}
		
		if(nrow(comment.pairs)>0)
		{
			for(pair in seq(nrow(comment.pairs)))
			{
				print(pair)
				rawtext = rawtext[-(comment.pairs[pair,1]:comment.pairs[pair,2])]
			}
		}
	}
	
	return (rawtext)
}

has.block <- function(finput,txt=NULL,blockname="characters2")
{
	retbool = (length(read.nexus.block(finput,txt,block=blockname,silent=T))!=0)
	return(retbool)	
}

has.characters2 <- function(finput,txt=NULL)
{
	return(has.block(finput,txt,"characters2"))
}



# Alternative way to read in characters
# Also works as a work around to readNexus' annoying habit
# or crashing when using type="data" but where trees block
# contains simmap formatted trees
read.characters2 <- function(finput,blockname="characters2")
{
	
	if(!file.exists(finput))
			stop("Assuming finput is a file and could not find it")
	
	rawtext = scan(finput,what=character(0),strip.white=T,sep="\n")		
	
	tmpfile = tempfile()
	tmphead = "#NEXUS\n"
	tmptaxa = c("BEGIN TAXA;",read.nexus.block(txt=rawtext,block="taxa"),"END;")
	tmptext = read.nexus.block(txt=rawtext,block=blockname)
	tmptext = c("BEGIN CHARACTERS;",tmptext,"END;")
	writeLines(c(tmphead,tmptaxa,tmptext),con=tmpfile)
	data2.part = readNexus(tmpfile,type="data",levels.uniform=F)
	rownames(data2.part) <- checkLabel(rownames(data2.part))
	return(data2.part)
}



