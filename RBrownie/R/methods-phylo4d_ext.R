
#-------------------------------------------------------
#  Methods 
#  -Are concerned with subnode stuff only.  
#   Higher-level tree modifications should be done
#  	with phylobase functions.
#-------------------------------------------------------



#--------------------------
# SIMMAP Processing Methods		
#--------------------------

# internal:
# get node names from newick string
get.nodenames<-function(newick.txt)
{
	nnames = character(0)
	ttmp = newick.txt
	nname.pat="(\\(|,|\\))([a-zA-Z0-9']{1,})(:|;)"
	junk = regexpr(nname.pat,ttmp)
	count = 0
	while(junk[1] != -1)
	{
		tmpname = substr(ttmp,(junk[1]+1),(junk[1]+attr(junk,"match.length")-2))
		nnames = append(nnames,tmpname)
		ttmp = sub(nname.pat,"",ttmp)	
		junk = regexpr(nname.pat,ttmp)
		count = count + 1
		stopifnot(count < 100000)
	}
	return(nnames)
}



# Check for evidence that this file contains simmap-formatted trees
#
is.simmap <- function(finput="",text=NULL)
{
	if(!is.null(text)){
		# TODO: check text for newlines and split on them if they exist.
		rawtext=text
	} else {
		if(!file.exists(finput))
			stop("Assuming finput is a file and could not find it")
		
		rawtext = scan(finput,what=character(0),strip.white=T,sep="\n")		
	}

	## TODO: split individual strings on ';' character
	
	# if there is only one string, then just check that
	if(length(rawtext)==1)
	{
		treesblock = rawtext
		treelines = 1
	} else {
		treesblock = read.nexus.block(txt=rawtext,block="trees",rm.comments=T)
		treelines = which(tolower(substr(treesblock,1,4))=="tree")
	}
	
	if(length(treesblock)==0)
	{
		warning("No trees in this file...")
		return(FALSE)
	}
	
	potentialsimmap = logical(length(treelines))
	count = 1
	for(linenumb in treelines)
	{
		# check for simmap style
		# remove first comment (should others be removed?
		junk =  gsub("\\[(.*?)\\]","",treesblock[linenumb])
		potentialsimmap[count] = as.logical(length(grep(":\\{.*?\\}",junk)))
		count = count + 1
	}
	
	return(any(potentialsimmap))
}



# read modified newick file
# citation: 
# Bollback J. P. (2006) SIMMAP: Stochastic character mapping of 
# discrete traits on phylogenies. BMC Bioinformatics. 7:88
# 
# Assume that branch lengths are present (otherwise, why use SIMMAP?)
# Assume that the root node can only have one simmap state (or, only use the first)
#
read.simmap <- function(file="",text=NULL, vers=1.1, ...)
{
	if(is.null(text))
		stop("Need to have text for now")

	# clear whitespace	
	text = gsub("\\s","",text)
	
	# add semicolon to end
	if(substring(text,nchar(text))!=";")
		text = paste(text,";",sep="")
	
	if(TRUE)
	{
		# add root node and internal node names
		count = 1
		while(regexpr("):",text)[1] != -1)
		{
			text = sub( "):", paste(")Internal",sprintf("%0.7d",count),":",sep=""), text)
			count = count + 1
		}
			
		# Root
		if(regexpr(");",text)[1] != -1)
			text = sub( ");", ")Root:0;" , text)
		
	}
	
	## Poor replacement for regular expressions
	#
	edge.ind = integer(0)
	edge.state = character(0)
	edge.len = numeric(0)
	junk = strsplit(text,"")[[1]]
	sub.branches = cbind( which(junk=="{"), which(junk=="}") )
	br.lens = numeric(nrow(sub.branches))
	
	for(ii in seq(nrow(sub.branches)))
	{
		br.total = 0
		# get the internal part
		within = paste(junk[seq( sub.branches[ii,1]+1 , sub.branches[ii,2]-1 )],collapse="")
		splitchar = ifelse(vers==1.0,";",":")
		within.sub = strsplit(within,splitchar)[[1]]
		for(jj in within.sub)
		{
			slptmp = strsplit(jj,",")[[1]]
			edge.ind = append(edge.ind,ii)
			edge.state = append(edge.state,slptmp[1])
			edge.len = append(edge.len,slptmp[2])
			br.total = br.total + as.numeric(slptmp[2])
		}
		br.lens[ii] = br.total
	}

	
	# horrible way to put it back together:
	# switch to regular expressions
	#
	newick.str = paste(junk,collapse="")
	for(ii in seq(nrow(sub.branches)))
	{
		replen = diff(sub.branches[ii,])+1
		substr(newick.str, sub.branches[ii,1], sub.branches[ii,2]) <- paste(rep(" ",replen),collapse="")
		substr(newick.str, sub.branches[ii,1], sub.branches[ii,2]) <- as.character(br.lens[ii])
	}
	newick.str = gsub("\\s","",newick.str) # strip whitespace
	
	# convert to ape format:
	tr = read.tree(text=newick.str)
	
	# make singleton nodes
	#
	edge.len = as.numeric(edge.len)
	all.labels = c(tr$tip.label,tr$node.label)
	
	# all names except the root node (which ought to be named, otherwise this will break)
	#
	#nnames = head(get.nodenames(newick.str),-1)  # This should match up with edge.ind
	nnames = get.nodenames(newick.str)
	scount = 1
	nTips = length(tr$tip.label)
	nInts = length(tr$node.label)
	
	# node information
	dataVal = character(0)
	dataNode = integer(0)
		
	# loop through split branches
	for(kk in unique(edge.ind))
	{
		# this information shouldn't change:
		is.root = F
		edgeind = which(edge.ind == kk)
		trind = which(nnames[kk] == all.labels)
		esplice = which(tr$edge[,2] == trind) # which edge to chop
		anc = tr$edge[esplice,1]
		dec = tr$edge[esplice,2]
	
		# assume that this is the root (no way to check for it otherwise..)
		if( length(esplice) == 0 )
		{
			dec = trind
			is.root = T
		}
		
		# store data information:
		dataNode = append(dataNode,dec)
		dataVal = append(dataVal, edge.state[edgeind[1]])

		newsingles = length(edgeind)-1
		
		if(newsingles==0 || is.root)
			next
		
		# add new singleton nodes
		for(mm in seq(newsingles))
		{
			#cat("- adding a new single\n")
			### add new internal nodes to 'phylo' object:
			#
			tr$Nnode = tr$Nnode + 1; nInts = nInts + 1
			tr$node.label = append(tr$node.label,sprintf("Singleton%0.7d",scount))
			nodeid = nTips + nInts # should be the last one available
			tr$edge[esplice,1] = nodeid
			tr$edge = rbind(tr$edge,c(anc,nodeid))
			tr$edge.length[esplice] = edge.len[edgeind[mm]]
			tr$edge.length = append(tr$edge.length,edge.len[edgeind[(mm+1)]])
			
			# store data information:
			dataNode = append(dataNode,nodeid)
			dataVal = append(dataVal, edge.state[edgeind[(mm+1)]])
			
			
			# update info:
			esplice = length(tr$edge.length) # should be the last one
			dec = nodeid
			scount = scount + 1
		}
		
	}
	
	# create phylo4d object
	new.labels = c(tr$tip.label,tr$node.label)
	tmpdf = data.frame("simmap_state"=dataVal, row.names=new.labels[dataNode])
	rettree = phylo4d(tr,all.data=tmpdf,missing.data="OK")
	
	#write.nexus(tr,file="written.tree")
	return(rettree)
}




# Read trees from a nexus file.  This function is only really necessary for 
# nexus files where trees have simmap formatting.  If they don't, then
# readNexus really should be used.
#
read.nexus.simmap <- function(finput="",text=NULL)
{
	if(!is.null(text)){
		# TODO: check text for newlines and split on them if they exist.
		rawtext=text
	} else {
		if(!file.exists(finput))
			stop("Assuming finput is a file and could not find it")
		
		rawtext = scan(finput,what=character(0),strip.white=T,sep="\n")		
	}

	treesblock = read.nexus.block(txt=rawtext,block="trees",rm.comments=T)
	
	if(length(treesblock)==0)
	{
		warning("No trees in this file...")
		return(FALSE)
	}

	## TODO: split individual strings on ';' character into newlines
	
	treelines = which(tolower(substr(treesblock,1,4))=="tree")
	if(length(treelines)==0)
		return(NULL)

	
	# make sure the the translated table is used if found
	transstart = which(sapply(treesblock,function(i) tolower(i)=="translate",USE.NAMES=F))
	translated = length(transstart)!=0
	transend=integer(0)
	taxatrans=""
	index=integer(0)
	taxname=""
	
	# assume that trees are not translated otherwise:
	# these should be in order numerically.
	if(translated)
	{
		transend = seq(transstart,min(treelines)-1)[min(grep(";",treesblock[seq(transstart,min(treelines)-1)]))]
		taxatrans = treesblock[(transstart+1):transend]
		
		# should ALWAYS be the last token:
		if(any(taxatrans==";"))
			taxatrans = taxatrans[-which(taxatrans==";")]
		
		index = as.integer(sapply(strsplit(taxatrans,"\\s"),function(i) i[1]))  # I'm assuming that they are integers
		taxname = sapply(strsplit(taxatrans,"\\s"),function(i) sub("(,|;)","",i[2]))
	}


	outtrees=list()
	count = 1
	for(linenumb in treelines)
	{
		# clearly, this will not work for non-simmap, non-newick trees
		tmpstr = tail(strsplit(treesblock[linenumb],"=")[[1]],1)
		
		# check for simmap style
		# remove first comment (should others be removed?
		if(is.simmap(text=treesblock[linenumb]))
		{
			trtmp = read.simmap(text=tmpstr)
		} else {
			trtmp = read.tree(text=tmpstr)
		}
		trtmp = as(trtmp,'phylo4')
		
		if(!inherits(trtmp,"phylo4")){
			cat("Trouble parsing line for trees:\n",treesblock[linenumb],"\n")
			stop()
		}
		
		outtrees = append(outtrees,trtmp)
		count = count + 1
	}
	
	# if the trees have been translated, then
	# add the names back in.
	#
	if(translated)
	{
		for(treeind in seq(length(outtrees)))
		{
			tnames = as.integer(tipLabels(outtrees[[treeind]]))
			if(!any(is.na(tnames)))
			{
				tipLabels(outtrees[[treeind]]) <- taxname[tnames]
			} else {
				warning("Tree translations found, but could not be mapped")
			}
		}
	}
	
	return(outtrees)
}




# expand singleton nodes into bifurcating nodes with one junk node
# this function is mainly for plotting and saving files
# TODO: convert all explicit calls to @'slot' to their accessor counterpart (e.g. tree@edge.length goes to edgeLength(tree)
#
expand.singles <- function(tree)
{
	# note: tips should be indexed 1...N, where N is the number of tips
	if(!is(tree,"phylo4"))
		stop("tree needs to be of class phylo4")
	
	if(hasSingle(tree))
	{
		tmptable=table(tree@edge[,1])
		snodes = as.integer(names(tmptable)[which(tmptable==1)])
		snodes = snodes[snodes!=0] # 0 is the root node (check on this...)
		
		nold = nrow(tree@edge)
		nnew = length(snodes)
		count = 1
		for(ii in snodes)
		{
			tree@label = append(sprintf("JUNK%0.7d",count),tree@label)
			tree@edge = rbind(tree@edge, c(ii, count))
			tree@edge.length = append(tree@edge.length,0)
			count = count + 1
		}
		
		# rearrange
		#tree@order = "unknown"
		rootind = which(tree@edge[,1] == 0)
		tree@edge[seq(nold),] = tree@edge[seq(nold),] + nnew
		tree@edge[seq(from=nold+1,to=nold+nnew),1] = tree@edge[seq(from=nold+1,to=nold+nnew),1] + nnew
		tree@edge[rootind,1] = 0
		
		# rename
		nnodes = nTips(tree) + nNodes(tree)
		names(tree@label) <- as.character(seq(1,nnodes))
		names(tree@edge.length) <- apply(tree@edge,1,paste,collapse="-")
	}
	
	return(tree)
}


# 
# sister function of expand.singles -> converts internal nodes with
# zero-length branches into singletons (in phylo4 format).  Will only work if
# the zero-len branch is connected to a tip
#
collapse.to.singles <- function(tree,by.name=NULL)
{
	if(is(tree,"phylo"))
		tree = as(tree,"phylo4")
	
	if(!is(tree,"phylo4"))
		stop("tree argument needs to be of class phylo4")

	# don't use the root node
	if(any(edgeLength(tree, seq(1,(nNodes(tree)+nTips(tree)))[nodeType(tree)!='root']) == 0))
	{
		# Remove zero-length branches and their tips
		zinds = which(edgeLength(tree)==0)
		torem = edges(tree)[zinds,]
		torem = torem[edges(tree)[zinds,1]!=0,][,2] # don't include root
		zinds = zinds[edges(tree)[zinds,1]!=0] 		# don't include root
		torem = torem[nodeType(tree)[torem]=="tip"]
		
		if(length(torem)!=0)
		{
			# remove excess tips
			rm.count = length(torem)
			tip.count = length(tipLabels(tree))
			total.count = length(labels(tree))
			
			tree@edge <- edges(tree)[-zinds,]  # TODO: add 'edges<-' to phylobase
			tree@edge.length = tree@edge.length[-zinds]
			tree@label = tree@label[-torem]
			
			## begin reindexing nodes
			# tips
			t.oldseq = seq(1, (tip.count))[-torem]
			t.newseq = seq(1,(tip.count-rm.count))
			t.replace.inds = which(tree@edge[,2] %in% t.oldseq)
			stopifnot(length(t.replace.inds) == length(t.oldseq))  # it should be a 1-1 relationship for tips
			tree@edge[,2] = replace(tree@edge[,2],t.replace.inds ,t.newseq)
			
			# internal
			n.oldseq = seq(tip.count+1,total.count)
			n.newseq = seq(tip.count+1-rm.count, total.count-rm.count)
			for(kk in seq(length(n.oldseq)))
			{
				replaceinds = which(tree@edge[,1] == n.oldseq[kk])
				if(length(replaceinds)!=0)
					tree@edge[replaceinds,1] = rep(n.newseq[kk],length(replaceinds))
					
				replaceinds = which(tree@edge[,2] == n.oldseq[kk])
				if(length(replaceinds)!=0)
					tree@edge[replaceinds,2] = rep(n.newseq[kk],length(replaceinds))
			}
			
			
			## rename
			names(tree@label) <- as.character(seq(1,length(labels(tree))))
			names(tree@edge.length) <- apply(tree@edge,1,paste,collapse="-")
		}
	} else {
		warning("No zero-length branches found to be removed")
	}
	return(tree)
}


# collapse singleton nodes using ape functions:
#
collapse.singletons <- function(phy)
{
	rettree = phy
	
	# if singletons exist
	if(hasSingle(rettree))
	{
		if(is(rettree,"phylo4d"))
		{
			tab=table(edges(rettree)[,1])
			snodeid = as.integer(names(tab)[which(tab==1)])
			snodeid = snodeid[snodeid!=0]  # 0 is a dummy node
			
			#
			ancs =  sapply(snodeid,function(i) which(edges(rettree)[,1] == i)) # where it is the ancestor
			decs =  sapply(snodeid,function(i) which(edges(rettree)[,2] == i)) # where it is the decendant
			newdata = data.frame(tdata(rettree,"all")[-decs,],row.names=labels(rettree)[-decs])
			colnames(newdata) <- colnames(tdata(rettree))
			
			# hack it for now:
			rettree = as(rettree,"phylo") # convert to ape
			rettree <- collapse.singles(rettree) # collapse singles
			rettree = phylo4d(rettree,all.data=newdata) # create new phylo4d object
		} else {
			rettree = as(rettree,"phylo") # convert to ape
			rettree <- collapse.singles(rettree) # collapse singles
			rettree = phylo4(rettree)
		}
	}
	
	return(rettree)
}






#------------------------------------------
# Modifying phylo4d-extension objects
# 
#------------------------------------------

# Internal:
# NOTE: labels are case-sensitive
#
.edge.index <- function(tree,anc,dec)
{
	if(!is(tree,'phylo4'))
		stop("Tree needs to inherit from class phylo4")
	
	# convert from characters to indices
	if(is.character(c(anc,dec))){
		anc = which(labels(tree) == anc)
		dec = which(labels(tree) == dec)
	}
			
	if(is.na(anc) || is.na(dec) || is.null(anc) || is.null(dec) || length(anc)!=1 || length(dec)!=1)
		stop("Must specify valid ancestor and decendent: ",anc,"-",dec)
	
	eind = which(edges(tree)[,1] == anc & edges(tree)[,2] == dec)
	if(length(eind) != 1){
		warning("No connection between ",anc," and ",dec,"\n")
		return(NA)
	}
	
	return(eind)
}


nSubNodes <- function(x)
{
	return (length(x@subnode.id))
}

hasSubNodes <- function(x)
{
	return (nSubNodes(x)!=0)
}

getSubNodeData <- function(x,colname)
{
	if(missing(colname))
		return(x@subnode.data)
	
	return(x@subnode.data[colname])
}

getSubNodePosish <- function(x)
{
	return(x@subnode.pos)
}


getSubNodeEdgeInds <- function(x)
{
	edgeinds = integer(nSubNodes(x))
	if(length(edgeinds)!=0)
	{
		for(ii in seq(nrow(x@subnode.branch)))
			edgeinds[ii] = .edge.index(x,x@subnode.branch[ii,1],x@subnode.branch[ii,2])
	}
	return (edgeinds)
}


# return empty data.frame
getEmptyDataFrame <- function(x)
{
	tmpdf = data.frame(x@data[0,])
	colnames(tmpdf) <- colnames(tdata(x))
	return(tmpdf)
}


# TODO: overload tdata here


# TODO: make this generic
# add a subnode
addSubNode <- function(x,anc,dec,position,dataf)
{
	eind = .edge.index(x,anc,dec)
	if(is.na(eind))
		stop("Failure to find edge from ",anc, " to ",dec,"\n")
	
	elen = edgeLength(x)[eind]
	if(elen == 0)
		stop("Cannot place subnode on a zero length branch")
	
	if(position > elen)
		stop("Position: ",position,", is greater than branch length: ",elen,"\n")
	
	# also check for overlapping subnodes
	
	# construct data frame:
	newdf = getEmptyDataFrame(x)
	if(is.data.frame(dataf))
	{
		if( all(names(dataf) == names(newdf)) )
			newdf = rbind(newdf, dataf)
		
	} else {

		if(ncol(newdf) == length(dataf)){
			newdf[1,] <- dataf
		
		} else {
			newdf[1,] <- rep(NA,ncol(newdf))
			# try to match up names
			ndf = names(dataf)
			count = 1
			for(nm in ndf){
				nmind = which(names(newdf)==nm)
				if(!is.null(nmind) && length(nmind) != 0)
					newdf[1,nmind] = dataf[count]
				count = count + 1
			}
		}
	}
	
	x@subnode.id = append(x@subnode.id, as.integer(nTips(x) + nNodes(x) + nSubNodes(x) + 1))
	x@subnode.data = rbind(x@subnode.data, newdf)
	x@subnode.branch = rbind(x@subnode.branch, edges(x)[eind,])
	x@subnode.pos = rbind(x@subnode.pos, rep(position,2))
	
	return (x)
}


#
showSubNodes <- function(x,ids)
{
	charlen = 80
	overlapchar = "*"
	# use them all by default:
	if(missing(ids))
		ids = x@subnode.id
	
	for(id in ids)
	{
		snid = which(x@subnode.id == id)
		anc = x@subnode.branch[snid,1]
		dec = x@subnode.branch[snid,2]
		tmpstr = rep("-",charlen)
		tmpstr[1] = 'o'; tmpstr[charlen] = 'o'
		eind = .edge.index(x,anc,dec)
		elen = edgeLength(x)[eind]
		breaksize = elen / charlen
		snpos = x@subnode.pos[snid,] * elen # Conrad: added this to get absolute position of subnode along the branch
		
		nbreaks = max(1, floor(diff(snpos) / breaksize) )
		from = floor(snpos[1] / breaksize)
		tmpstr[seq(from,length.out=nbreaks)] = rep(overlapchar,nbreaks)
		tmpstr = paste(tmpstr,collapse="")
		names(tmpstr) <- sprintf("Subnode at ~%0.2f on branch: %d to %d (brlen=%0.2f)",snpos[1]+diff(snpos)/2,anc,dec,elen)
		print(tmpstr)
		cat("\n")
	}
}


#--------------------------
# Plotting
# TODO: virtually all of it.
#--------------------------
plotPhyext <- function(x,y,...)
{
	# TODO: figure out a way to place node data over nodes...
	cat("Plotting extensions\n")
	treePlot(x,plot.data=F, ...)
}

setMethod('plot', signature(x='phylo4d_ext', y='missing'), function(x, y, ...) {
    plotPhyext(x, ...)
})



#-----------------------------
# Set Generics
#-----------------------------

# show
setMethod("show","phylo4d_ext", function(object){ printphylo4(object); showSubNodes(object)})

## snData
setGeneric("snData", function(x, ...) {
    standardGeneric("snData")
})

## snData<-
setGeneric("snData<-", function(x, ..., value) {
    standardGeneric("snData<-")
})


setMethod("snData", signature(x="phylo4d_ext"),
  function(x) {
	return(x@subnode.data)
})


