#---------------------------------------------
# Process Nexus files  
# -	Extra methods for extracting different information
# 	from nexus-formatted files.
#
#---------------------------------------------


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


# get tree weights from file or tree string
# Assuming this format: 
# [.... &W -122235 ........]
#
get.tree.weights <- function(finput,text=NULL,splitchar=" ")
{
	if(!is.null(text)){
		# TODO: check text for newlines and split on them if they exist.
		rawtext=text
	} else {
		if(!file.exists(finput))
			stop("Assuming finput is a file and could not find it")
		
		rawtext = scan(finput,what=character(0),strip.white=T,sep="\n")		
		rawtext = read.nexus.block(txt=rawtext,block="trees")
	}
	
	comments = get.nexus.comments(text=rawtext)
	tokens = strsplit(comments,splitchar)
	find.weight <- function(ii,ww = "&W")
	{
		return(ii[(which(toupper(ii) == ww))+1])
	}
	weighs = unlist(lapply(tokens,find.weight))
	
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
					newtokens[curline] = paste(tmpline,char,sep="")
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
	
	start.ind = agrep(paste("begin",tolower(blockname)),filetext,ignore.case=T)
	if(length(start.ind) == 0)
		return (integer(0))
	
	end.ind = grep("end;",filetext,ignore.case=T)
	end.ind = end.ind[head(which(end.ind > start.ind),1)]

	return( c((start.ind),(end.ind)) )
}



# Internal function to parse / check assumptions block
.process.assumptions <- function(obj,block.txt)
{
	
	if(!is(obj[[1]],"brownie"))
		stop("Processed object needs to be of class brownie")
	
	
	for(aline in block.txt)
	{
		tokens = strsplit(aline,"\\s")[[1]]
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
			nm = paste("TAXSET",taxname,sep="_")
			taxaI = data.frame(all=rep(0,nTips(obj[[1]])))
			names(taxaI) <- nm
			
			for(tind in seq(length(obj)))
			{
				# convert if needed
				if(!inherits(obj[[tind]],"phylo4d"))
					obj[[tind]] = phylo4d(obj[[tind]])
				
				taxaI[,1] = sapply(tipLabels(obj[[tind]]),function(i) ifelse(i %in% taxinds,1,0),simplify=T)
				#names(taxaI) <- nm
				obj[[tind]] = addData(obj[[tind]],tip.data=taxaI)
				
			}
		}
	}
	
	return (obj)
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
		datatypes = rep(genericData(),ndatcols)
		datatypes[sapply(seq(ndatcols),is.numeric)] = contData()
		datatypes[sapply(seq(ndatcols),is.factor)] = discData()
		datatypes[grep("TAXSET_",names(datvals))] = taxaData()
		
		for(tind in seq(length(obj)))
		{
			obj[[tind]]@datatypes <- datatypes
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


# Internal:
# Get text of a read nexus block
read.nexus.block<-function(finput,txt=NULL,block)
{	
	if(!is.null(txt)){
		# Using the text argument is not recommended
		rawtext=txt
	} else {
		
		if(!file.exists(finput))
			stop("Assuming finput is a file and could not find it")
		
		rawtext = scan(finput,what=character(0),strip.white=T,sep="\n")		
	}
	
	inds = .get.nexus.block.inds(blockname=block,text=rawtext)
	if(length(inds)==0)
	{
		warning(paste("This file has",block, "no block"))
		return (character(0))
	}
	
	# TODO: split up newlines if they exist
	return (rawtext[(inds[1]+1):(inds[2]-1)])
}



