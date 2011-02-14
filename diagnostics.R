

# Method which compare two trees for equality
#
cmptrees <- function(t1,t2,extended=TRUE)
{
	retval = TRUE
	
	if(nNodes(t1) != nNodes(t2)) {
		warning("nNodes don't match")
		retval = FALSE
	}
	
	if(nTips(t1) != nTips(t2)) {
		warning("nTips don't match")
		retval = FALSE
	}
	
	if(nEdges(t1) != nEdges(t2)) {
		warning("nEdges don't match")
		retval = FALSE
	}
	
	if(xor(hasSingle(t2),hasSingle(t1))) {
		warning("hasSingle don't match")
		retval = FALSE
	}
	
	if(xor(isRooted(t2),isRooted(t1))) {
		warning("isRooted don't match")
		retval = FALSE
	}
	
	# Not sure if this will always be case:	
	if(rootNode(t1) != rootNode(t2)) {
		warning("rootNode don't match")
		retval = FALSE
	}
	
	# edge lengths:
	if(any(is.na(edgeLength(t1))))
		t1@edge.length[is.na(edgeLength(t1))] <- 0.000
		
	if(any(is.na(edgeLength(t2))))
		t2@edge.length[is.na(edgeLength(t2))] <- 0.000
	
	if(any(sort(round(edgeLength(t1),2)) != sort(round(edgeLength(t2),2)))) {
		warning("rootNode don't match")
		retval = FALSE
	}	
	
	if(extended)
	{
		# compare subnodes:
		retval = retval && cmpSubNodes(t1,t2)
		retval = retval && cmpData(t1,t2)
	}
	retval
}


cmpSubNodes<-function(t1,t2)
{
	retval=TRUE
		
	if(xor(hasSubNodes(t2),hasSubNodes(t1))) {
		warning("hasSubNodes don't match")
		retval = FALSE
	}

	# If they both have subnodes:
	if(hasSubNodes(t1) && hasSubNodes(t2))
	{
		for(ii in seq(ncol(sndata(t1))))
		{
			if(!is.factor(sndata(t1)[,ii]))
			{
				if(!all(sort(sndata(t1)[,ii]) == sort(sndata(t2)[,ii]))) {
					warning("sndata don't match",ii)
					retval = FALSE
				}
			} else {
				#commonlevs = intersect(levels(sndata(t1)[,ii]),levels(sndata(t2)[,ii]))
				if(!all(sort(as.character(sndata(t1)[,ii])) == sort(as.character(sndata(t2)[,ii])))) {
					warning("sndata don't match",ii)
					retval = FALSE
				}
			}
		}
	
		if(!all(round(sort(apply(snposition(t1),1,mean))) == round(sort(apply(snposition(t2),1,mean))))) {
			warning("snposition don't match")
			retval = FALSE
		}

		# This is not necessarily the same between trees (node ids could be different for similar trees)
		#
		#if(!all(sort(apply(snbranch(t1),1,sum)) == sort(apply(snbranch(t2),1,sum)))) {
		#	warning("snbranch don't match")
		#	retval = FALSE
		#}

		if(!all(sort(edgeLength(t1)[getSubNodeEdgeInds(t1)])==sort(edgeLength(t2)[getSubNodeEdgeInds(t2)])))
		{
			warning("snbranch don't match")
			retval = FALSE			
		}
	
	}
	retval
}


cmpData <- function(t1,t2)
{
	retval = TRUE

	if(xor(hasData(t2),hasData(t1))) {
		warning("hasData does not match")
		retval = FALSE
	}
	
	if(hasData(t1) && hasData(t2))
	{
		datanames = colnames(tdata(t1))
		
		for(ii in seq(ncol(tdata(t1))))
		{
			if(!is.factor(tdata(t1)[,ii]))
			{
				if(!all(sort(tdata(t1)[,ii]) == sort(tdata(t2)[,ii]))) {
					warning("tdata don't match",ii)
					retval = FALSE
				}
			} else {
				
				tmp1 = tdata(t1)[,ii]
				tmp1 = tmp1[which(!is.na(tmp1))]
				tmp1 = as.character(tmp1)
				tmp1 = tmp1[which(tmp1!="NA")]

				tmp2 = tdata(t2)[,ii]
				tmp2 = tmp2[which(!is.na(tmp2))]
				tmp2 = as.character(tmp2)
				tmp2 = tmp2[which(tmp2!="NA")]

				if( !all(sort(tmp1) == sort(tmp2)) ) {
					warning("tdata don't match (factor)",ii)
					retval = FALSE
				}
			}
		}
	}
	
	return(retval)
}

# compare the output of two data.frames from Brownie analyses:
#
cmpAnalysis<-function(df1,df2,allow.similar=TRUE,dferror=0.05)
{
	retval = TRUE
	df1names = colnames(df1)
	df2names = colnames(df2)
	df1mods = df1$Model
	df2mods = df2$Model
	
	for(ii in seq(length(df1mods)))
	{
		df1ind = which(df1mods == df1mods[ii])
		df2ind = which(df2mods == df2mods[ii])

		for(jj in seq(length(df1[df1ind,])))
		{
			colind1 = jj
			colind2 = which(df2names == df1names[colind1])

			# compare values:
			if(is.na(df1[df1ind,colind1]) && is.na(df2[df2ind,colind2]))
			{	
				# do nothing
			} else {
			
				is.same = (df1[df1ind,colind1] == df2[df2ind,colind2])

				if(!is.same && allow.similar && !is.character(df1[df1ind,colind1]))
				{
					val1 = df1[df1ind,colind1]
					val2 = df2[df2ind,colind2]
					upper.val = val1 + (val1*dferror)
					lower.val = val1 - (val1*dferror)
					if( !(val2 <= upper.val && val2 >= lower.val) )
					{
						warning("Values did not match for column '", df1names[colind1],"' (",colind1," ",colind2,") on rows ",df1ind," ",df2ind ) 
						retval = FALSE
					}
				} else {
					if(!is.same && !allow.similar){
						warning("2 Values did not match for column '", df1names[colind1],"' on row ",df1ind) 
						retval = FALSE				
					}
				}
					
			}
		}
	}

	retval
}





