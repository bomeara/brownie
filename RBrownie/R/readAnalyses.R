#--------------------------------------------------
# Run and return results of Phylogenetic analyses 
# 
#--------------------------------------------------

# Method to scan a brownie output file and return it as a dataframe
# Brownie generally output multiple header lines, these need to be 
# removed before the dataframe is returned.
#
read.analysis.output <- function(filename,txt=NULL,rowsep='\n',colsep='\t')
{
	rettab = character(0)
	if(!is.null(txt))
	{
		rettab = txt
	} else {
		# use the file:			
		ff = filename
		if(!file.exists(ff))
			stop("Could not find the file specified")
		
		rettab = scan(ff,what="character",sep=rowsep,strip.white=T)
	}
	
	headercol=1
	headers = strsplit(rettab[headercol],colsep)[[1]]
	datacols= integer(0) # put out all data columns
	count = 1
	for (line in rettab)
	{
		if(line != rettab[headercol]){
			datacols = c(datacols,count)
		}
		count = count + 1
	}
	dfout = data.frame(matrix(NA,nrow=length(datacols),ncol=length(headers)))
	colnames(dfout) <- headers
	datasep = strsplit(rettab[datacols],'\t')
	for(jj in seq(length(datasep)))
	{
		dfout[jj,] = datasep[[jj]] 
	}
	
	# post conditioning (changing to numeric):
	oldwarn = options()$warn; options(warn=-1) # suppress these warnings
	for(kk in seq(ncol(dfout)))
	{
		if(!any(is.na(as.numeric(dfout[,kk]))))
			dfout[,kk] = as.numeric(dfout[,kk])
	}
	options(warn=oldwarn)
	
	return(dfout)
	
}

# 
# ellipsis stuff if passed to read.analysis.output
# 
read.discrete.output <- function(filename,txt=NULL,...)
{
	warnpattern="^WARNING:"
	output=character(0)
	if(!is.null(txt))
	{
		output = txt
	} else {
		output = scan(filename,what="character",sep="\n",strip.white=T)
	}
	
	warnlines = grep(warnpattern,output)
	if(length(warnlines)!=0)
	{
		cat("Analysis returned some warnings:\n")
		print(output[warnlines])
		output = output[-warnlines]
	}
	
	ret = read.analysis.output(txt=output)
	return(ret)
}



