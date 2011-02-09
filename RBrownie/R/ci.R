# get coalescent intervals from tip dated trees


tipdate.ci<-function(tr,show.plot=F)
{
	
	if ( !(class(tr) %in% c("multiPhylo","phylo")) )
	  	stop("object \"tr\" is not of class \"phylo\"")
	if(!is.binary.tree(tr))
		stop("Need a binary tree")
	if(!is.rooted(tr))
		stop("Need a rooted tree")

	nleaves = length(tr$tip.label)
	ileaves = seq(1,nleaves)
	#rnode = nleaves + 1 # root node
	rnode =  which(tabulate(tr$edge[,2])==0)[1]
	cipos<<-numeric(0)
	citype<<-numeric(0)
	cinode<<-numeric(0)
	
	isleaf<-function(n)
	{
		return( (length(which(tr$edge[,1]==n)) == 0) )
	}
	
	addlen<-function(len)
	{
		cipos<<-c(cipos,len)
	}
	
	addtype<-function(type)
	{
		citype<<-c(citype,type)
	}	
	
	addnode<-function(nnode)
	{
		cinode<<-c(cinode,nnode)
	}
	
	# USE SAPPLY for recursion from now on!
	# Start at root node and trace back from there....
	# node = node number,
	# tr = tree (should stay constant)
	# len = interval lengths
	#
	chknode<-function(node,len)
	{
		#Sys.sleep(0.25)
		#cat("node: ",node," length: ",len,"\n")
		# check for leaf node
		if(isleaf(node))
		{
			addlen(len)
			addtype(0)		
			addnode(node)
		} else {
			addlen(len)
			addtype(1)	
			addnode(node)
			# tr$edge 	 col1 = FROM,
			#		 col2 = TO
			inds = which(tr$edge[,1]==node)
			if(show.plot) edgelabels(as.character(tr$edge.length[inds]),inds,adj=-0.5)			
			Recall(tr$edge[inds[1],2],(len + tr$edge.length[inds[1]]))
			Recall(tr$edge[inds[2],2],(len + tr$edge.length[inds[2]]))
		}
	}
	
	if(show.plot) plot.phylo(tr,show.tip.label=F)

	chknode(rnode,0)
	sorder = order(cipos, decreasing=TRUE)
	tmppos = cipos[sorder]
	tmptype = citype[sorder]
	#rm(list=c("cipos","citype","cinode"))
	
	lin = numeric(0) #lineages
	il = abs(diff(tmppos)) #interval.length
	ic = length(cipos)-1 # interval.count
	td = diff(range(cipos)) #total.depth
	
	count=0;
	for (i in tmptype)
	{
		if(i==0) count = count + 1 	# sampling event
		else count= count - 1		# coalescent event
		lin = append(lin,count)
	}
	
	return(list(lineages=lin[1:(length(lin)-1)],interval.length=il,interval.count=ic,total.depth=td,I=tmptype[-1]))
}


setGeneric("treeHeight", function(x) { standardGeneric("treeHeight") })
setMethod("treeHeight", signature(x="phylo4d_ext"),
  function(x) {	
	options(warn=-1)  # suppress warnings for this sections
	retval = tipdate.ci(as(x,"phylo"))$total.depth
	options(warn=0)  
	return(retval)
})

setMethod("treeHeight", signature(x="phylo4d"),
  function(x) {
	options(warn=-1)  # suppress warnings for this sections
	retval = tipdate.ci(as(x,"phylo"))$total.depth
	options(warn=0)  
	return(retval)
})

setMethod("treeHeight", signature(x="phylo4"),
  function(x) {
	options(warn=-1)  # suppress warnings for this sections
	retval = tipdate.ci(as(x,"phylo"))$total.depth
	options(warn=0)  
	return(retval)
})

