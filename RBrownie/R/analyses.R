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



