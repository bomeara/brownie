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
	
	if(length(rettab)==0)
	{
		warning("Return table has no characters.  Looks like brownie did not return anything or execute possibly.")
		return(data.frame(NULL))
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
	if(length(datasep)!=0)
	{
		for(jj in seq(length(datasep)))
		{
			dfout[jj,] = datasep[[jj]] 
		}
	} else {
		warning("No information could be retrieved from analyis. Brownie might not have run properly.")
		return(dfout)
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

# Read continuous test output
read.continuous.output <- function(filename,txt=NULL,...)
{
	# The first line with tabs is the header:
	tabpattern="\t"
	output=character(0)
	if(!is.null(txt))
	{
		output = txt
	} else {
		output = scan(filename,what="character",sep="\n",strip.white=T)
	}
	
	tablines = grep(tabpattern,output)
	if(length(tablines)!=0)
	{
		output = output[tablines]
	} else {
		stop("Failed to find any output that could be coersed into a table\nOUTPUT:\n",output)
	}
	
	ret = read.analysis.output(txt=output)
	
	return(ret)
}


# Internal function:
scan.textout <- function(output,strip.white=T)
{
	# scan text out and remove empty lines:
	strsplit(output,"\n")->tmp
	tmp = unlist(tmp)
	emptylines = which(nchar(tmp)==0)
	if(length(emptylines)!=0)
		tmp = tmp[-emptylines]
	
	# remove whitespace
	if(strip.white)
	{
		tmp = sub("\\s+$","",tmp) # trim trailing whitespace
		tmp = sub("^\\s+","",tmp) # trim leading whitespace
	}
	
	return(tmp)
}


# Internal: read treesout into a list of phylogenetic tree objects
# -this process will not vary between analyses
#
scan.treesout <- function(output)
{
	tmp = list()
	for(ii in seq(length(output)))
	{
		tmpstr = gsub("'","",output[ii])  # brownie puts little tip marks as internal node names
		if(is.simmap(text=tmpstr))
		{
			tmptree = read.simmap(text=tmpstr)
		} else {
			tmptree = read.tree(text=tmpstr)
		}
		tmp = append(tmp,phyext(tmptree))
	}
	return(tmp)
}



#-----------------------------------------------------
# Summary functions (for each brownie output)
# 
#-----------------------------------------------------

summary.ratetest <- function(ratedf,txt=NULL,short=FALSE)
{
	cat("Summary of ratetest results:")
	bootstrapped=FALSE
	# convert text into data.frame
	if(!is.null(txt))
	{
		ratedf = read.continuous.output(txt=scan.textout(txt))
	}

	headers = names(ratedf)
	if(is.null(ratedf$Tree) || is.null(ratedf$Char))
		stop("Could not find Trees or Char columns in the ratedf dataframe.\nAvailable columns are:",headers)
		
	# check for output columns named 'param', these
	# are only(?) set when reps!=0 or, to put it another way,
	# when bootstrapping is done.
	#
	if(length(grep("^param",headers,value=T)) != 0)
		bootstrapped = TRUE
	
	utrees = unique(ratedf$Tree)
	uchars = unique(ratedf$Char)
	tcgrid = expand.grid(utrees,uchars)
	
	if(!short)
	{
		for(ii in seq(nrow(tcgrid)))
		{
			rowind = which(ratedf$Tree == tcgrid[ii,1] & ratedf$Char == tcgrid[ii,2])
			if(length(rowind)!=1)
				stop("Could not find unique row index for tree,char combo (",tcgrid[ii,],")")
			
			cat("\n====================================================\n\n")
			cat("Tree =",tcgrid[ii,1],", character =",tcgrid[ii,2],"\n")
			cat("----------------------------\n\n")
			cat("Results for the inverse of the taxaset specified (NOT):\n")
			cat(sprintf("Anc state = %14.6f\nRate = %14.6f\n-lnL = %14.6f\n\n",ratedf$anc_1[ii],ratedf$rate_1[ii],ratedf$' -lnL_1'[ii]))
			cat("Results for the taxaset specified:\n")
			cat(sprintf("Anc state = %14.6f\nRate = %14.6f\n-lnL = %14.6f\n\n",ratedf$anc_2[ii],ratedf$rate_2[ii],ratedf$' -lnL_2'[ii]))
			cat("----------------------------\n\n")
			cat("Single rate support (model A) vs. Multiple rate support (model B)\n\n")
			cat(sprintf("            -lnL          AIC          AICc                 \n"))
			cat(sprintf("A: %14.6f %14.6f %14.6f\n",ratedf$" -lnL_A"[ii],ratedf$AIC_A[ii],ratedf$AICc_A[ii]))
			cat(sprintf("B: %14.6f %14.6f %14.6f\n",ratedf$" -lnL_B"[ii],ratedf$AIC_B[ii],ratedf$AICc_B[ii]))
			cat(sprintf("diff:             %14.6f %14.6f\n",ratedf$"AIC dif"[ii],ratedf$"AICc dif"[ii]))
			cat("Chi-squared p-value:",ratedf$"chi p"[ii],"\n")
			if(bootstrapped)
				cat("Bootstrap p-value:",ratedf$"param p"[ii],"\n")
			cat("----------------------------\n")
		}
	}
	cat("\nOverall summary of results:\n")
	cat("(B/b = strong/weak support for multiple rate models)\n")
	cat(sprintf("%s\t%s\t%s\t%s\t%s\t%s\n","tree","char","AIC","AICc","Chi",ifelse(bootstrapped,"param","")))
	for(ii in seq(nrow(tcgrid)))
	{
		rowind = which(ratedf$Tree == tcgrid[ii,1] & ratedf$Char == tcgrid[ii,2])
		# assume it's the last column:
		vals = ratedf[,ncol(ratedf)][ii]
		vals = strsplit(vals,"")[[1]]
		cat(sprintf("%d\t%d\t%s\t%s\t%s\t%s\n",ratedf$Tree[ii],ratedf$Char[ii],vals[1],vals[2],vals[3],ifelse(bootstrapped,vals[4],"")))
	}

}




summary.discrete <- function()
{
	
}




summary.cont <- function()
{
	
}


