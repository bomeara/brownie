#-----------------------
# Phylo4d class extended
# @ author Conrad Stack
#
# This extension will support "sub-nodes" so that branches can be fragmented to 
# accomadate the results of ancestral state reconstructions (state changes along 
# branches).
#-----------------------

# @param subnode.id  	ID for this node 
# @param subnode.data  	data.frame the should look like superclass data.frame
# @param subnode.branch The branch that this subnode is located on. (Duplicate information - the reason this uses node id's instead of edge indices is due to the fact that edge indices seem unstable)
# @param subnode.pos  	The distance to the subnode from the ancestor (root) of the branch. Can be a non-zero range
setClass("phylo4d_ext", 
			representation(	subnode.id="integer",
							subnode.data="data.frame",
							subnode.branch="matrix",  
							subnode.pos="matrix",
							weight="numeric" ),
			,contains="phylo4d",
			prototype=prototype(
				subnode.id=integer(0),
				subnode.data=data.frame(NULL),
				subnode.branch=matrix(0,nrow=0,ncol=2),
				subnode.pos=matrix(0,nrow=0,ncol=2),
				weight=numeric(0)
				)
			)


#------------------
# check validity
#------------------
validPhylo4d_ext <-function(object)
{
	
	retval = TRUE

	# number of ids and branches are the same
	retval = retval && (length(object@subnode.id) == nrow(object@subnode.branch))

	# number of branches and number of entries in the data matrix are the same
	if(length(object@subnode.data)!=0)
		retval = retval && (nrow(object@subnode.branch) == nrow(object@subnode.data))
	
	# number of ids and number of rows in position matrix are similar
	retval = retval && (length(object@subnode.id) == nrow(object@subnode.pos))
	
	# Check that branch ids exist in superclass 
	if(nrow(object@subnode.branch)!=0)
		for(br.ind in seq(nrow(object@subnode.branch)))
			retval = retval && any(apply(edges(object),1,function(i) all(object@subnode.branch[br.ind,] == i)))
		
	# check that data.frames in 'data' slot and 'subnode data' slot are structurally similar
	# retval = retval && all( names(object@data) == names(object@subnode.data) )
	
	return (retval)
}
setValidity("phylo4d_ext",validPhylo4d_ext)


#----------------------------------------------
#  Constructors: extended phylo4d  
#  (overloaded)
#----------------------------------------------

# make it a generic function with return class of phylo4d_ext (or list of these)
setGeneric("phyext", function(x, ...) { standardGeneric("phyext")}, valueClass=c("phylo4d_ext","list") )

# return itself (otherwise it will use the function below):
setMethod("phyext", "phylo4d_ext",
    function(x, ...) {
	    return(x)
})
	
# if first arg is a phylo4d object:
setMethod("phyext", "phylo4d",
    function(x, snode.data=NULL, snode.branch=NULL, snode.pos=NULL, ...) {
	
	# default values:
	if(missing(snode.data) || is.null(snode.data))
		snode.data = data.frame(x@data[0,])
	
	if(missing(snode.branch) || is.null(snode.branch))
		snode.branch = matrix(0,nrow=0,ncol=2)
	
	if(missing(snode.pos) || is.null(snode.pos))
		snode.pos = matrix(0,nrow=0,ncol=2)
	
	# convert to matrix is it isn't already
	if(is.numeric(snode.pos))
		snode.pos = matrix(rep(snode.pos,2),ncol=2)
	
	if(is.numeric(snode.branch))
		snode.branch = matrix(snode.branch,ncol=2)
		
	# Convert any singletons into subnodes:
	#
	if(hasSingle(x))
	{
		
		###########
		# Changed (7/11) by conrad to process more than one subnode:
		#
		# process subnodes:		
		
		
		# Get singleton node ids:
		tab=table(edges(x)[,1])
		snodeid = as.integer(names(tab)[which(tab==1)])
		snodeid = snodeid[snodeid!=0]  # 0 is a dummy node
		
		
		# Get singleton indices: (edges where singletons are the descendant)
		#sings = sapply(snodeid,function(i) which(edges(x)[,2] == i))
		sings = snodeid
		
		# Get singleton data: (why are edge indices used here(sings) instead of node ids?)
		#decs.data = data.frame(tdata(x,"all")[sings,])
		#dnames <- colnames(tdata(x))
		#names(decs.data) <- dnames
		decs.data = tdata(x,"all")[sings,,drop=F]
		
		
		########################################################################################
		# Setup other temporary data containers for transferring singleton data to subnodes
		# 
		# Idea:
		# 1.)	einds gets a first singletons:
		# 		new.ancestor <<----------->> first.singleton
		# 2.)	For loop the processes all subsequent singletons, looping until non-singleton is found
		# 3.) 	Put data back together
		#
		
		# 1.)
		einds = which(apply(edges(x),1,function(i) !(i[1]%in%snodeid) && i[2]%in%snodeid  )) # get all first singletons-on-a-branch
		alle = edges(x)[einds,]
		if(!is.matrix(alle))  alle = matrix(alle,nrow=1)
		elens = edgeLength(x)[einds]
		#edata = data.frame()
		#colnames(edata) <- dnames
		edata = tdata(x,"all")[alle[,2],,drop=F]
		ee = cbind(alle,elens)
		newnodes = matrix(0,ncol=3,nrow=0)
		
		# 2.)
		for(ii in seq(nrow(ee)))
		{
			anc = ee[ii,1]
			dec = ee[ii,2]
			
			# if decendent is still a singleton:
			if(dec %in% snodeid)
			{
				newlens = ee[ii,3]
				newdata = data.frame(NULL)
	
				# travel back until a non-singleton is found:
				while(dec %in% snodeid)
				{
					#cat("dec=",dec," - ")
					
					# length to new subnode:
					tmplen = unname(edgeLength(x)[which(edges(x)[,1]==dec)])
					newlens = append(newlens, (tail(newlens,1)+tmplen) )
					#cat(levels(tmpdata[1,])[as.integer(tmpdata[1,])],"\n")
					
					# new subnode data:
					#tmpdata = data.frame(decs.data[which(sings==dec),])
					#colnames(tmpdata) <- dnames
					tmpdata = decs.data[which(sings==dec),,drop=FALSE]
					newdata = rbind( newdata, tmpdata )						
					
					dec = edges(x)[which(edges(x)[,1]==dec),2]
				}
				
				# use any
				newlens = head(tail(newlens,-1),-1)
				newdata = tail(newdata,-1)
				if(length(newlens)!=0)
				{
					for(kk in seq(length(newlens)))
					{
						newnodes = rbind(newnodes,c(anc,dec,newlens[kk]))
					}
					edata = rbind(edata,newdata)
				}
				ee[ii,2] <- dec
			}
		}
		
		# 3.)
		eetmp = rbind(ee,newnodes)
		ancs = unname(eetmp[,1])
		decs = unname(eetmp[,2])
		snode.length = unname(eetmp[,3])
		snode.data = edata
		#names(snode.data) <- dnames
		
		
		###############################################################################################################
		
		# new labels:
		ancs.label = labels(x)[ancs]
		decs.label = labels(x)[decs]
		
		# Remove singletons from original tree
		x = collapse.singletons(x)
		
		# condition the branches:
		for(ii in seq(length(ancs)))
		{
			# NOTE: This function assumes that labels are unique (might not always be true...)
			snode.branch = rbind(snode.branch, edges(x)[.edge.index(x, ancs.label[ii], decs.label[ii]),])
			# Conrad: make subnode position a fraction of the parent branch length:
			#snode.pos = rbind(snode.pos, rep(snode.length[ii],2))
			snode.pos = rbind(snode.pos, rep(snode.length[ii],2)/edgeLength(x)[.edge.index(x,snode.branch[ii,1],snode.branch[ii,2])]) 
		}
	}
	
	# create dummy ids for subnodes 
	# NOTE: subnode.ids not currently being used
	#
	idvals = integer(0)
	if(nrow(snode.branch)!=0)
	{
		idstart = (nTips(x) + nNodes(x) + 1)
		idvals = as.integer(seq(from = idstart, length.out=nrow(snode.branch)))
	} 
	
	# Annotate collapsed tree with subnode data:
	retobj = new("phylo4d_ext",
		x,
		subnode.id=idvals, 
		subnode.data = snode.data,
		subnode.branch = snode.branch, 
		subnode.pos = snode.pos)
	
	return(retobj)
})



setMethod("phyext", "phylo4",
    function(x, snode.data=NULL, snode.branch=NULL, snode.pos=NULL, ...) {

	# create new phylo4 object (allow arguments to be passed via ellipsis)
	phyd = phylo4d(x,...)
	
	phyext(phyd,snode.data,snode.branch,snode.pos)
})



setMethod("phyext", "phylo",
    function(x, snode.data=NULL, snode.branch=NULL, snode.pos=NULL, ...) {
	phyd = phylo4d(x,...)
	phyext(phyd,snode.data,snode.branch,snode.pos)
	
})


setMethod("phyext","list",
	function(x,...){
	
	x = sapply(x,phyext,...)
	return(x)
})


# assume character points to a file or is a piece of text
# TODO: check if trees are simmap formatted!
setMethod("phyext", "character",
    function(x, snode.data=NULL, snode.branch=NULL, snode.pos=NULL, ...) {
	 
	# Figure out what x's format is:
	tr = NULL
	if(file.exists(x))
	{
		# x points to a file  
		if(is.simmap(finput=x)){
			tr = try(read.nexus.simmap(x)[[1]],silent=T)
		} else {
			tr = try(read.tree(x),silent=T)
		}
		
		if(is(tr,"try-error"))
		{
			tr = try(read.nexus(x),silent=T)
			if(is(tr,"try-error"))
				stop("If first argument is a file then it must contains trees in newick or nexus format")
		}
	} else {
		# x is a character string:
		if(is.simmap(text=x)){
			tr = try(read.simmap(text=x),silent=T)
		} else {
			tr = try(read.tree(text=x),silent=T)
		}
		if(is(tr,"try-error"))
			stop("If first argument is a string, then it must be newick-formated.")
	}
	
	phyd = phylo4d(tr,...)
	phyext(phyd,snode.data,snode.branch,snode.pos)
})






