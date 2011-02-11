
#--------------------------------------------------
# Read / Write brownie 
# -	This file relies heavily on processNexus.R and 
#	methods-phylo4d_ext.R
#
#--------------------------------------------------

# readBrownie
# writeBrownie

setGeneric("writeBrownie", function(x,...) { standardGeneric("writeBrownie")} )


readBrownie<-function(fname)
{
	
	if(!file.exists(fname))
		stop(paste("File",fname,"cannot be found in",getwd()))
	
	filetxt = scan(fname,what=character(0),strip.white=T,sep="\n")
	brownie.part = read.nexus.block(txt=filetxt,block="BROWNIE",rm.comments=T)
	assumptions.part = read.nexus.block(txt=filetxt,block="ASSUMPTIONS",rm.comments=T)
	
	issim = is.simmap(fname)
	if(!issim){
		
		# this function should read in character data
		# wrap in list to make it compatible with read.nexus.simmap output
		phy.part = readNexus(fname,levels.uniform=F)
		if(!is.list(phy.part))
		{
			phy.part = list(phy.part)
		}
		
		# rename:
		for(jj in seq(length(phy.part)))
			tipLabels(phy.part[[jj]])<-checkLabel(tipLabels(phy.part[[jj]]))
		
	} else {
		
		# should return a list of phylo4 objects with singleton nodes
		# NOTE (2/11) - this now returns 'phylo4d_ext' objects
		phy.part = read.nexus.simmap(fname) 
		
		# process data part (the hard way)
		#data.part = readNexus(fname,type="data")
		data.part = read.characters2(fname,blockname="characters")  # work around; added 08/04
		
		if(length(data.part) > 0)
		{
			# Different versions of phylobase work differently:
			warning("Having to use sketchy regular expressions to add data (i.e. The difficult way)")
			data.names = tolower(rownames(data.part))
			
			for(tind in seq(length(phy.part)))
			{
				# NOTE: doing this reordering is necessary because phylobase does not seem to
				# 		return unadulterated taxa names where reading character data.
				# convert to phylo4d 
				if(!inherits(phy.part[[tind]],"phylo4d"))
					phy.part[[tind]] = phylo4d(phy.part[[tind]])
				
# 				if(any(data.names != tolower(tipLabels(phy.part[[tind]]))))
# 				{
# 					
				#tipmods = tolower(sub("[^[:alnum:]]","",tipLabels(phy.part[[tind]]),extended=T))
				tipmods = tolower(sub("[^[:alnum:]]","",tipLabels(phy.part[[tind]])))
				data.names = tolower(sub("[^[:alnum:]]","",data.names))
				
				if(any(sort(tipmods)!=sort(data.names)))
					stop("Could not match up tree names against data names")
				
				neworder=unname(sapply(tipmods,function(i) which(i == data.names)))
				data.part.tmp = data.part[neworder,,drop=F]
					
# 				} else {
# 					# This will probably never be called until NCL is updated
# 					data.part.tmp = data.part
# 				}

				phy.part[[tind]] = addData(phy.part[[tind]],tip.data=data.part.tmp,match.data=F)
				phy.part[[tind]] = phyext(phy.part[[tind]])  # NOTE (2/11): this should be redundant since read.nexus.simmap now returns only 'phylo4d_ext' objects
			}
		}
	}
	
	if(!inherits(phy.part[[1]],"phylo4"))
		stop("Object of class phylo4 was not created.  Email author of this shoddy code.")
	
	brau.new = list()
	for(tind in seq(length(phy.part)))
	{
		 brau.new = append(brau.new,new("brownie",phy.part[[tind]],commands=brownie.part))
	}

	brau.new = .process.datatypes(brau.new)  # may or may not be useful
	brau.new = .process.assumptions(brau.new,assumptions.part,issim)
	
	
	#
	# Read characters 2 if it exists:
	if(has.characters2(fname))
	{
		cat("Processing CHARACTERS2 block\n")
		data2.part = read.characters2(fname)
# 		for(xx in seq(length(brau.new)))
# 			brau.new[[xx]] = addData(brau.new[[xx]],data2.part)		
		
		brau.new <- addData(brau.new,data2.part)
	}
	
	# read in weights if they exist:
	# BUG: bad wieghts (negative values?) are crashing the program
	#if(has.weights(fname))
	if(FALSE)
	{
		cat("Processing tree weights...")
		weights.part = get.tree.weights(fname)
		if(length(weights.part)==0){
			cat("not found\n")
		}else{
			weight(brau.new) <- weights.part
			cat("done\n")
		}
	}
	
	if(length(brau.new)==1)
		return(brau.new[[1]])
	
	return(brau.new)
}


.write.brownie.block <- function(phytree)
{
	bblock=character(0)
	cmds = commands(phytree)
	cmds[grep(";$",cmds)] = paste(cmds[grep(";$",cmds)],";",sep="")
	if(length(cmds)!=0)
	{
		bblock = paste("begin brownie;",paste(commands(phytree),collapse="\n"),"end;",sep="\n")
	}
	bblock
}


.write.taxa.strings <- function(phytree)
{
	if(!is(phytree,"phylo4d"))
		stop("phytree needs to have a data section.")
	
	taxsetnames = colnames(taxasets(phytree))
	tstring = character(0)
	
	if(length(taxsetnames)!=0)
	{
		for(ii in seq(length(taxsetnames)))
		{
			tmpname = sub("^TAXSET_","",taxsetnames[ii])
			tmp = paste(taxa.charvect(phytree,taxsetnames[ii]),collapse=" ")
			tmpname = paste("taxset",tmpname)
			tstring = c(tstring, paste(paste(tmpname,tmp,sep="="),";",sep="") )
		}
	}
	tstring
}


# write assumptions
#
.write.assumptions.block <- function(phytree)
{
	assout = character(0)
	if(!is(phytree,"phylo4d"))
		stop("phytree needs to have a data section.")
	
	tstr = .write.taxa.strings(phytree)
	if(length(tstr)!=0)
	{
		tstr = paste(tstr,collapse="\n")
		assout = paste("BEGIN assumptions;",tstr,"END;",sep="\n")
	}
	
	return (assout)
}

# return conditioned character string
#
.write.characters.block <- function(xdf,blockname="CHARACTERS",dtype,missing.char="?")
{	
	# set up state labels:
	use.state.labels=F
	state.labels=character(0)
	
	# For disc data, figure out what the symbols should be 
	discover.symbols<-function(dat,zero.based=T)
	{			
		syms = integer(0)
		if(ncol(dat)!=0){
			for(ii in seq(ncol(dat)))
			{
				syms = c(syms,unique((as.integer(dat[,ii]))))
			}
		}
		syms = sort(unique(syms))
		
		# assume that syms are in order
		if(zero.based)
		{
			if(syms[1] != 0){
				offs = -syms[1]
				syms = syms + offs
			}
		}
		
		return(syms)
	}
	
	# convert factors to zero-based ints
	# CAUTION: this function might not work in next NCL version
	convert.to.int <- function(dat,zero.based=T)
	{
		offset = 0
		if(zero.based)
			offset = -1
		
		for(ii in seq(ncol(dat)))
		{
			dat[,ii] = as.integer(dat[,ii])+offset
		}
		return(dat)
	}
	
	# generate state labels:
	make.state.labels <- function(dat)
	{
		cnames = colnames(dat)
		outlabs = character(ncol(dat))
		for(ii in seq(ncol(dat)))
		{
			# order by
			snames = unique(as.character(dat[,ii])[order(as.integer(dat[,ii]))])
			outlabs[ii] = sprintf("%d %s / %s",ii,cnames[ii],paste(snames,collapse=" "))
		}
		
		return(outlabs)
	}
	
	
	if(!is.data.frame(xdf))
		stop("Internal function .write.characters.block needs a data.frame as the first argument")
	
	header = paste("BEGIN ",blockname,";",sep="")
	header.title = paste("TITLE ",blockname,"_matrix;",sep="")
	header.dims = sprintf("DIMENSIONS NTAX=%d NCHAR=%d;",nrow(xdf),ncol(xdf))
	header.format = sprintf("FORMAT DATATYPE=%s MISSING=%s",ifelse(dtype==contData(),"CONTINUOUS","STANDARD"),missing.char)   # TODO: add GAP, SYMBOLS 
	if(dtype == discData()){
		# This check is any of the levels are NOT integers
		# if they are integers, then it is assumed that they do not 
		# need to be writen as state labels
		#
		use.state.labels = any(is.na(as.integer(levels(xdf[,1]))))
		if(use.state.labels){
			state.labels = make.state.labels(xdf)
			xdf = convert.to.int(xdf) # this does the zero-based conversion
		}
		header.format = sprintf("%s SYMBOLS=\"%s\";",header.format,paste(discover.symbols(xdf),collapse=" "))
	} else {
		header.format = paste(header.format,";",sep="")
	}
	
	if(use.state.labels){
		header.labels = sprintf("CHARSTATELABELS\n\t%s;", paste(state.labels,collapse=","))
	}else{
		header.labels = sprintf("CHARSTATELABELS\n\t%s;", paste(paste(seq(ncol(xdf)),colnames(xdf)),collapse=","))
	}

	mmatrix = "MATRIX"
	#mmatrix.data = unname(cbind(rownames(xdf),apply(xdf,2,as.character)))
	mmatrix.data = apply(xdf,2,as.character)
	if(any(is.na(mmatrix.data)))
		mmatrix.data[which(is.na(mmatrix.data),T)] <- missing.char
	
	mmatrix.data = apply(mmatrix.data,1,paste,collapse=" ")
	mmatrix.data = unname(cbind(rownames(xdf),mmatrix.data))
	mmatrix.data = unname(apply(mmatrix.data,1,paste,collapse="\t"))
	mmatrix.end = ";\n\nEND;"
	

	return(c(header,
			header.title,
			header.dims,
			header.format,
			header.labels,
			mmatrix,
			mmatrix.data,
			mmatrix.end))

}


setMethod("writeBrownie",signature(x="brownie"),
	function(x,file=NULL,rmsimmap=TRUE) {
		return( writeBrownie(list(x),file=file,rmsimmap=rmsimmap) )
})

# write nexus file with trees and characters
# TODO: -Use text streams instead of temp files
#		-Better way to convert CR/LF (in Windows)
#
setMethod("writeBrownie", signature(x="list"),
	function(x, file=NULL, rmsimmap=TRUE) {
		
		# temporary files for nexus blocks:
		#
		tmp1 = tempfile()  # TREES / TAXA
		tmp2 = tempfile()  # CHARACTERS
		tmp3 = tempfile()  # CHARACTERS2
		tmp4 = tempfile()  # BROWNIE 
		tmp5 = tempfile()  # ASSUMPTIONS
		
		success = TRUE
		datablock = FALSE
		datablock2 = FALSE
		brownieblock = FALSE
		assumptionsblock=FALSE
		
		#if(missing(usechar) || is.null(usechar))
		#	usechar = names(tdata(x,"tip"))[1]
		
		# Perpare tree
		#phy = as(x[[1]],'phylo')
		write.nexus.simmap(x,file=tmp1)
		
		# if there is tip data to be written
		#
		if(hasTipData(x[[1]]))
		{
			if(rmsimmap)
			{
				x[[1]] = rmdata(x[[1]],'simmap_state')
			}
			
			dtypes = datatypes(x[[1]])
			if(any(dtypes == genericData()))
				warning("Excluding undefined datatypes.")
		
			udtypes = sort(unique(dtypes[which(dtypes %in% c(contData(),discData()))]),T)
			
			if(length(udtypes) > 0)
			{
				
				# Perpare first CHARACTERS block:
				#
				dat = tdata(x[[1]],"tip")[,which(dtypes==udtypes[1]),drop=F]
				tnames = rownames(dat)	
				#write.nexus.data(dat,file=tmp2,"STANDARD",datablock=F)
				tmpss = .write.characters.block(dat,"CHARACTERS",udtypes[1])
				writeLines(tmpss,tmp2)
				datablock=TRUE
				
				if(length(udtypes)==2)
				{
					dat = tdata(x[[1]],"tip")[,which(dtypes==udtypes[2]),drop=F]
					#write.nexus.data(dat,file=tmp3,"STANDARD",datablock=F)
					tmpss = .write.characters.block(dat,"CHARACTERS2",udtypes[2])
					writeLines(tmpss,tmp3)
					datablock2=TRUE
				}
			}
		}
				
		
		if(hasCommands(x[[1]])){
			tmpss = .write.brownie.block(x[[1]])
			if(length(tmpss)==0) stop("Brownie commands found, but could not be written.  Email authors")
			writeLines(tmpss,tmp4)
			brownieblock = TRUE
		}
		
		if(hasTaxasets(x[[1]]))
		{
			tmpss = .write.assumptions.block(x[[1]])
			if(length(tmpss)==0) stop("Taxa sets found, but could not be written.  Email authors")
			writeLines(tmpss,tmp5)
			assumptionsblock=TRUE
		}
		
		# write to specified file:
	
		success = file.copy(tmp1,file,overwrite=T)
		
		if(datablock)
			success = success && file.append(file,tmp2)
		if(datablock2)
			success = success && file.append(file,tmp3)
		if(assumptionsblock)
			success = success && file.append(file,tmp5)
		if(brownieblock)
			success = success && file.append(file,tmp4)
		
		if(!success)
			stop("Failed to copy between temporary files")
		
		# IF windows, convert to windows path:
		# (brownie crashes otherwise)
		#
# 		if(.Platform$OS.type=="windows")
# 		{
# 			tmpconvert = tempfile()
# 			file = gsub("\\\\","/",file)
# 			tmpconvert = gsub("\\\\","/",tmpconvert )
# 			sysstr = paste("tr -d '\\015' < ", file, " > ", tmpconvert)
# 			shell(sysstr)
# 			success = success && file.copy(tmpconvert,file,overwrite=TRUE)
# 			if(!success)
# 				stop("Failed to convert line endings on this Windows machine")
# 			
# 			unlink(tmpconvert)
# 		} 

		# alternative
		if(.Platform$OS.type=="windows")
		{
			file = gsub("\\\\","/",file)
			strs = readLines(file)
			writeLines(strs,con=file,sep="\r")
		} 
		
		
		# delete temporary files:			
		unlink(tmp1)
		unlink(tmp2)
		unlink(tmp3)
		unlink(tmp4)
		unlink(tmp5)
		
		# return
		return( success )
})
	
