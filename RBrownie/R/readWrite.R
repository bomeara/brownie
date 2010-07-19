
#--------------------------------------------------
# Read / Write brownie 
# -	This file relies heavily on processNexus.R and 
#	methods-phylo4d_ext.R
#
#--------------------------------------------------


readBrownie<-function(fname)
{	
	if(!file.exists(fname))
		stop(paste("File",fname,"cannot be found in",getwd()))
	
	filetxt = scan(fname,what=character(0),strip.white=T,sep="\n")
	brownie.part = read.nexus.block(txt=filetxt,block="BROWNIE")
	assumptions.part = read.nexus.block(txt=filetxt,block="ASSUMPTIONS")
	
	if(!is.simmap(fname)){
		
		# this function should read in character data
		# wrap in list to make it compatible with read.nexus.simmap output
		phy.part = list(readNexus(fname))  
		
	} else {
		
		# should return a list of phylo4 objects with singleton nodes
		phy.part = read.nexus.simmap(fname) 
		
		# process data part (the hard way)
		data.part = readNexus(fname,type="data")
		
		if(length(data.part) > 0)
		{
			warning("Having to use sketchy regular expressions to add data (i.e. The difficult way)")
			data.names = tolower(rownames(data.part))
			
			# TODO: fix phylobase!
			for(tind in seq(length(phy.part)))
			{
				# NOTE: doing this reordering is necessary because phylobase does not seem to
				# 		return unadulterated taxa names where reading character data.
				# convert to phylo4d 
				if(!inherits(phy.part[[tind]],"phylo4d"))
					phy.part[[tind]] = phylo4d(phy.part[[tind]])
				
				tipmods = tolower(sub("[^[:alnum:]]","",tipLabels(phy.part[[tind]]),extended=T))
				neworder=unname(sapply(tipmods,function(i) which(i == data.names)))
				data.part.tmp = data.part[neworder,]
				phy.part[[tind]] = addData(phy.part[[tind]],tip.data=data.part.tmp,match.data=F)
				phy.part[[tind]] = phyext(phy.part[[tind]])
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
	brau.new = .process.assumptions(brau.new,assumptions.part) 
	brau.new = .process.datatypes(brau.new)  # may or may not be useful
	
	if(length(brau.new)==1)
		return(brau.new[[1]])
	
	return(brau.new)
}


.write.brownie.block <- function(phytree)
{
	outfile = tempfile()
	cat(paste("begin brownie;",paste(phytree@commands,collapse="\n"),"end;",sep="\n"),file=outfile)
	return(outfile)
}

# write nexus file with trees and characters
# NOTE: this converts file to Unix line-endings
#		need to figure out a better way to do this. 
#
write.nexus.both <- function(phytree,file="",usechar=NULL)
{
	brownieblock=FALSE
	retbool = TRUE
	
	if(is(phytree,'phylo4d'))
	{	
		if(is(phytree,'brownie'))
			brownieblock = TRUE
		
		if(missing(usechar) || is.null(usechar))
			usechar = names(tdata(phytree,"tip"))[1]
			
		phy = as(phytree,'phylo')
		dat = tdata(phytree,"tip")[usechar]
		dnames = rownames(dat)
		dat = as.character(dat[,1])
		names(dat) <- dnames
		
		tmp1 = tempfile()
		tmp2 = tempfile()
		tmp4 = ""
		
		write.nexus(phy,file=tmp1)
		write.nexus.data(dat,file=tmp2,"STANDARD",datablock=F)
		retbool = retbool && file.append(tmp1,tmp2)
		
		if(brownieblock){
			tmp4 = .write.brownie.block(phytree)
			retbool = retbool && file.append(tmp1,tmp4)
		}
		
		if(.Platform$OS.type=="windows")
		{
			tmp3 = tempfile()
			# IF windows:
			# convert to windows path:
			tmp1 = gsub("\\\\","/",tmp1)
			tmp3 = gsub("\\\\","/",tmp3)
			sysstr = paste("tr -d '\\015' < ", tmp1, " > ", tmp3)
			shell(sysstr)
			retbool = retbool && file.copy(tmp3,file,overwrite=TRUE)
			unlink(tmp3)
		} else {
			retbool = retbool && file.copy(tmp1,file,overwrite=TRUE)
		}
		
		unlink(tmp1)
		unlink(tmp2)
		
		if(brownieblock)
			unlink(tmp4)
		
	} else {
		warning("Object phytree did not inherit from phylo4d")
	}
	
	return( retbool )
}

