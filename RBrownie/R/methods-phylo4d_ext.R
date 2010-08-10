
#-------------------------------------------------------
#  Methods 
#  -Are concerned with subnode stuff only.  
#   Higher-level tree modifications should be done
#  	with phylobase functions.
#-------------------------------------------------------

# Table of Contents:

## SIMMAP Methods:
# get.nodenames
# is.simmap
# read.simmap
# read.nexus.simmap
# expand.singles
# collapse.to.singles
# collapse.singletons
# write.simmap

## Subnode Methods:
# .edge.index
# nSubNodes
# hasSubNodes
# getSubNodeData
# getSubNodePosish
# getSubNodeEdgeInds
# getSubNodeEdgeInds
# getEmptyDataFrame
# addSubNode
# showSubNodes

## Generics
setGeneric("sndata", function(x, ...) { standardGeneric("sndata") })
setGeneric("sndata<-", function(x,datnames=NULL, value) { standardGeneric("sndata<-") })
setGeneric("snid", function(x, ...) { standardGeneric("snid") })
setGeneric("snposition", function(x, ...) { standardGeneric("snposition") })
setGeneric("snbranch", function(x, ...) { standardGeneric("snbranch") })
setGeneric("rmdata", function(x,index) { standardGeneric("rmdata")} )
setGeneric("weight", function(x) { standardGeneric("weight")} )
setGeneric("weight<-", function(x,value) { standardGeneric("weight<-")} )
setGeneric("hasWeight",function(x,strict=TRUE) { standardGeneric("hasWeight")} )

# subnode
setGeneric("hasSubNodes", function(x) { standardGeneric("hasSubNodes") })

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
	outtrees = NULL
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

	count = 0
	for(linenumb in treelines)
	{
		print(linenumb)
		# clearly, this will not work for non-simmap, non-newick trees
		tmpstr = tail(strsplit(treesblock[linenumb],"=")[[1]],1)
		
		# check for simmap style
		# remove first comment (should others be removed?
		if(is.simmap(text=treesblock[linenumb]))
		{
			trtmp = unname(read.simmap(text=tmpstr))
		} else {
			trtmp = read.tree(text=tmpstr)
		}
		if(!is(trtmp,'phylo4'))
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



# write modified newick file:
# -This is kind of a hack: it basically writes subnodes and lengths to a new character label and then
#  creates a new 'phylo' class using only those labels (without edge lengths).  Then it uses
#  APEs algorithm to write those labels to a newick string.  
#
# -If any of the SIMMAP datatypes are is.na, then they are left out!  TODO:  Look into changing this in the future
#
write.simmap <- function(x,usestate="simmap_state",file="",vers=1.1,...)
{
	splitchar = ifelse(vers==1.0,";",":")
	if(hasSubNodes(x))
	{
		tdat = tdata(x)[,usestate,drop=F]
		sdat = sndata(x)[,usestate,drop=F]
		snedges.inds = apply(snbranch(x),1,function(i) .edge.index(x,i[1],i[2]))
		es = edges(x)[,2]
		elens = edgeLength(x)
		newlenlab=character(length(es))
		for(ii in seq(length(es)))
		{
			nodeid=es[ii]
			if(!(ii %in% snedges.inds))
			{
				if(!is.na(tdat[nodeid,1]))
					newlenlab[ii] = paste(tdat[nodeid,1] ,",", elens[ii],sep="")
			} else {
				snind = which(snedges.inds == ii)
				snlens = snposition(x)[snind,] * elens[ii]
				snstates = sdat[snind,1]
				if(!is.matrix(snlens))
					snlens = matrix(snlens,nrow=1)
				snpos = apply(snlens,1,mean)
				snpos = sort(snpos,T)
				newlenlab[ii] = paste(tdat[nodeid,1],",",(elens[ii]-max(snpos)),sep="")
				
				if(length(snpos)>1)
					for(jj in seq(length(snpos)-1))
						newlenlab[ii] = paste(newlenlab[ii],splitchar,snstates[jj],",",(snpos[jj]-snpos[(jj+1)]),sep="")
				
				newlenlab[ii] = paste(newlenlab[ii],splitchar,tail(snstates,1),",",tail(snpos,1),sep="")
				
			}
		}
		oldlabs = labels(x)[es]
		names(newlenlab) <- oldlabs
		oldlabs[which(is.na(oldlabs))] <- ""
		newlab = paste(oldlabs,":{", newlenlab ,"}",sep="")
		newlab[which(newlenlab=="")] <- ""  # remove any <NA> data
		
		# reorder:
		newlab = newlab[order(es)]
		ntype = nodeType(x)
		phy = as(x,'phylo')
		newedges = edges(x)
		newphy = list(edge=newedges,tip.label=newlab[which(ntype=="tip")],node.label=newlab[which(ntype!="tip")],Nnode=nrow(newedges))
		class(newphy) <- "phylo"
		
		
		################################################################
		# borrowed code from APE (write.tree.R):
		#
		output.tree.names=FALSE
		append = FALSE
		digits = 10
		brl <- !is.null(newphy$edge.length)
		nodelab <- !is.null(newphy$node.label)
		f.d <- paste("%.", digits, "g", sep = "")
		cp <- function(s) STRING <<- paste(STRING, s, sep = "")
		add.internal <- function(i) {
		    cp("(")
		    br <- which(newphy$edge[, 1] == i)
		    for (j in br) {
		        desc <- newphy$edge[j, 2]
		        if (desc > n) add.internal(desc)
		        else add.terminal(j)
		        if (j != br[length(br)])  cp(",")
		    }
		    cp(")")
		    if (nodelab) cp(newphy$node.label[i - n])
		    if (brl) {
		        cp(":")
		        cp(sprintf(f.d, newphy$edge.length[which(newphy$edge[, 2] == i)]))
		    }
		}
		add.terminal <- function(i) {
		    cp(newphy$tip.label[newphy$edge[i, 2]])
		    if (brl) {
		        cp(":")
		        cp(sprintf(f.d, newphy$edge.length[i]))
		    }
		}
		n <- length(newphy$tip.label)
		STRING <- if (output.tree.names) paste("ERROR!", "(", sep = "") else "("
		br <- which(newphy$edge[, 1] == n + 1)
		for (j in br) {
		    desc <- newphy$edge[j, 2]
		    if (desc > n) add.internal(desc)
		    else add.terminal(j)
		    if (j != br[length(br)]) cp(",")
		}
		if (is.null(newphy$root.edge)) {
		    cp(")")
		    if (nodelab) cp(newphy$node.label[1])
		    cp(";")
		} else {
		    cp(")")
		    if (nodelab) cp(newphy$node.label[1])
		    cp(":")
		    cp(sprintf(f.d, newphy$root.edge))
		    cp(";")
		}
		if (file == "") return(STRING)
	    cat(STRING, file = file, append = append, sep = "\n")
	    
		################################################################		
	} else {
		phy = as(x,'phylo')
		write.tree(phy,file,...)
	}
}	


# This is the main write function for phylo4d_ext
# Mainly ripped from APE write.nexus function
#
#
write.nexus.simmap <- function(obj, file = "", translate = TRUE)
{
	if(!is.list(obj))
	{
		if(!is(obj,"phylo4d_ext"))
			stop("This function is only made to work with phylo4d_ext objects.")
		
		obj <- list(obj)
	} else {
		if(!all(sapply(obj,is,'phylo4d_ext')))
			stop("This function is only made to work with phylo4d_ext objects or lists of such.")
	}
	ntree <- length(obj)
	
    cat("#NEXUS\n", file = file)
    cat(paste("[R-package APE, ", date(), "]\n\n", sep = ""),file = file, append = TRUE)
	
    #N <- length(obj[[1]]$tip.label)
	N <- length(tipLabels(obj[[1]]))

        cat("BEGIN TAXA;\n", file = file, append = TRUE)
        cat(paste("\tDIMENSIONS NTAX = ", N, ";\n", sep = ""),file = file, append = TRUE)
        cat("\tTAXLABELS\n", file = file, append = TRUE)
        cat(paste("\t\t", tipLabels(obj[[1]]), sep = ""),sep = "\n", file = file, append = TRUE)
        cat("\t;\n", file = file, append = TRUE)
        cat("END;\n", file = file, append = TRUE)
   	
    cat("BEGIN TREES;\n", file = file, append = TRUE)
    if (translate) {
        ## We take arbitrarily the labels of the first tree, and
        ## translate them as "1", "2", "3", ...
        cat("\tTRANSLATE\n", file = file, append = TRUE)
        tmp <- checkLabel(tipLabels(obj[[1]]))
        X <- paste("\t\t", 1:N, "\t", tmp, ",", sep = "")
        ## We remove the last comma:
        X[length(X)] <- gsub(",", "", X[length(X)])
        cat(X, file = file, append = TRUE, sep = "\n")
        cat("\t;\n", file = file, append = TRUE)
        token <- as.character(1:N)
        names(token) <- tipLabels(obj[[1]])
        tipLabels(obj[[1]]) <- token
        if (ntree > 1) {
            for (i in 2:ntree)
                tipLabels(obj[[i]]) <- token[tipLabels(obj[[i]])]
            class(obj) <- NULL
        }
    } else {
        for (i in 1:ntree)
          tipLabels(obj[[i]]) <- checkLabel(tipLabels(obj[[i]]))
    }
    for (i in 1:ntree) {
	    tprefix = "\tTREE * UNTITLED"
	    #weights
	    if(hasWeight(obj[[i]])){
		    tprefix = sprintf("%s [&W %f]",tprefix,weight(obj[[i]]))
	    }
	    #is rooted
        if (isRooted(obj[[i]])){
        	#cat("\tTREE * UNTITLED = [&R] ", file = file, append = TRUE)
        	tprefix = paste(tprefix," = [&R] ",sep="")
    	}else{
	    	#cat("\tTREE * UNTITLED = [&U] ", file = file, append = TRUE)
	    	tprefix = paste(tprefix," = [&U] ",sep="")
		}
		cat(tprefix, file = file, append = TRUE)
        cat(write.simmap(obj[[i]], file = ""),"\n", sep = "", file = file, append = TRUE)
    }
    cat("END;\n", file = file, append = TRUE)
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


setMethod("hasSubNodes", signature(x="phylo4d_ext"),
  function(x) {
	return(nSubNodes(x)!=0)
})


setMethod("hasSubNodes", signature(x="phylo4"),
  function(x) {
	return(FALSE)
})

setMethod("hasSubNodes", signature(x="phylo"),
  function(x) {
	return(FALSE)
})



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

# return empty data.frame styled like
# @data slot
getEmptyDataFrame <- function(x)
{
	tmpdf = data.frame(x@data[0,])
	colnames(tmpdf) <- colnames(tdata(x))
	return(tmpdf)
}


# TODO: overload tdata here

# 
addSubNode <- function(x,anc,dec,position,dataf)
{
	
	if(!is(x,'phylo4d_ext')){
		warning("x is not an extended phylo4d object")
		return(x)
	}
	
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
showSubNodes <- function(x)
{
	if(hasSubNodes(x))
	{
		charlen = 80
		brchar = "-"
		terminalchar = "-"
		regchar="0"
		overlapchar = "*"
		# use them all by default:
		
		einds = apply(x@subnode.branch,1,function(i) .edge.index(x,i[1],i[2]))
		
		for(eind in unique(einds))
		{
			snid = which(einds == eind)
			anc = x@subnode.branch[snid,1]
			dec = x@subnode.branch[snid,2]
			
			# setup output string
			tmpstr = rep(brchar,charlen)
			tmpstr[1] = terminalchar
			tmpstr[charlen] = terminalchar
			
			# get relative positions of nodes:
			elen = edgeLength(x)[eind]
			breaksize = elen / charlen
			snpos = x@subnode.pos[snid,] * elen 
			if(!is.matrix(snpos)) 
				snpos = matrix(snpos,nrow=1)
			
			for(xx in seq(length(snid)))
			{
				nbreaks = max(1, floor(diff(snpos[xx,]) / breaksize) )
				from = floor(snpos[1] / breaksize)
				tmpstr[seq(from,length.out=nbreaks)][tmpstr[seq(from,length.out=nbreaks)]==regchar] = overlapchar
				tmpstr[seq(from,length.out=nbreaks)][tmpstr[seq(from,length.out=nbreaks)]==brchar] = regchar
				#tmpstr[seq(from,length.out=nbreaks)] = rep(regchar,nbreaks)
			}
			
			# collapse and print:
			tmpstr = paste(tmpstr,collapse="")
			if(length(snid)==1){
				names(tmpstr) <- sprintf("Subnode at ~%0.2f on branch: %d to %d (brlen=%0.2f)",snpos[1]+diff(snpos[1,])/2,anc,dec,elen)
			} else {
				names(tmpstr) <- sprintf("%d Subnodes on branch: %d to %d (bren=%0.2f)",length(snid),anc[1],dec[1],elen)
			}
			print(tmpstr)
			cat("\n")
		}
		print("0 indicates a subnode; * indicates 2+ subnodes overlapping. Positions are relative.")
	}
}


#-----------------------------
# Set Generics
#-----------------------------

# show
setMethod("show","phylo4d_ext", function(object){ printphylo4(object); showSubNodes(object)})

setMethod("snid", signature(x="phylo4d_ext"),
  function(x) {
	return(x@subnode.id)
})

setMethod("snid", signature(x="list"),
	function(x) {
		return(x[[1]]@subnode.id)
})

setMethod("sndata", signature(x="phylo4d_ext"),
  function(x) {
	return(x@subnode.data)
})

setMethod("sndata", signature(x="list"),
  function(x) {
	return(x[[1]]@subnode.data)
})


# this adds data
setReplaceMethod("sndata", signature(x="phylo4d_ext"),
  function(x,datnames=NULL,value) {
	
	if(!is.data.frame(value) && (length(value) != nSubNodes(x)))
	{
	  warning("Can only add vectors of data if they are the same length as the number of subnodes")
	  return(x)
	}
	if(!is.data.frame(value)){
		value = data.frame(value)
	}
	if(!is.null(datnames))
	{
		if(is.character(datnames) && length(datnames)==ncol(value))
		{
			names(value) <- datnames
		} else {
			warning("Not using specified names")
		}
	}
	x@subnode.data = cbind(x@subnode.data,value)
	return(x)
})


setReplaceMethod("sndata",signature(x="list"),
	function(x,datnames=NULL,value) {
		for(ii in seq(length(x)))
			sndata(x[[ii]],datnames=datnames) <- value
			
	return(x)
})


setMethod("rmdata", signature(x="phylo4d_ext",index="numeric"),
  function(x,index) {
	  if(length(index)>0 && index <= ncol(tdata(x)))
		{
			x@data = x@data[,-index,drop=F]
			if(hasSubNodes(x))
				x@subnode.data = x@subnode.data[,-index,drop=F]
		}
	
	return(x)
})


setMethod("rmdata", signature(x="phylo4d",index="numeric"),
  function(x,index) {
	  if(length(index)>0 && index <= ncol(tdata(x)))
		x@data = x@data[,-index,drop=F]
	
	return(x)
})

setMethod("rmdata", signature(x="phylo4d",index="character"),
  function(x,index) {
	inds = which(names(tdata(x)) %in% index)
	return(rmdata(x,inds))
})


setMethod("rmdata", signature(x="list",index="ANY"),
	function(x,index) {
		x = sapply(x,rmdata,index)
	return(x)
})


setMethod("snposition", signature(x="phylo4d_ext"),
  function(x) {
	return(x@subnode.pos)
})

setMethod("snposition", signature(x="list"),
  function(x) {
	return(x[[1]]@subnode.pos)
})

setMethod("snbranch", signature(x="phylo4d_ext"),
  function(x) {
	return(x@subnode.branch)
})

setMethod("snbranch", signature(x="list"),
  function(x) {
	return(x[[1]]@subnode.branch)
})


#---------------
## WEIGHT:
#
#---------------
setMethod("weight", signature(x="phylo4d_ext"),
  function(x) {
	return(x@weight)
})

setMethod("weight", signature(x="list"),
	function(x) {
		if(hasWeight(x))
			return( sapply(x,weight) )
		
		return(numeric(0))
})

setReplaceMethod("weight", signature(x="phylo4d_ext"),
  function(x,value) {
	x@weight = value
	return(x)
})

setReplaceMethod("weight",signature(x="list"),
	function(x,value) {
	
	if(length(x) != length(value))
		stop("Replacement values need to be same length as the list")
	
	for(ii in seq(length(x)))
		weight(x[[ii]]) <- value[ii]
	
	return(x)
})

setMethod("hasWeight",signature(x="phylo4d_ext"),
	function(x,strict=TRUE){
		return( length(x@weight)!=0 )
})


setMethod("hasWeight",signature(x="list"),
	function(x,strict=TRUE){
		retbool=ifelse(strict,TRUE,FALSE)
		for(ii in seq(length(x))){
			tmpbool = hasWeight(x[[ii]])
			if(strict && !tmpbool)
				return(FALSE)
			
			if(!strict && tmpbool)
				return(TRUE)
			
		}
		return(retbool)
})

# 
# setMethod("addData", signature(x="phylo4d_ext"),
# 	function(x,...,snode.data=NULL) {
# 		
# 		oldcols = names(tdata(x))
# 		# Add data the normal way:
# 		x = getMethod("addData","phylo4d")(x,...)
# 		
# 		newcols = setdiff(names(tdata(x)),names(sndata(x)))
# 		supercols = setdiff(oldcols,names(tdata(x)))
# 		
# 		cat("newcols: ", newcols,"\n")
# 		
# 		# Add data to subedges:
# 		if(is.null(snode.data))
# 		{
# 			if(length(newcols)!=0)
# 			{
# 				newdat = data.frame(matrix(NA,nrow=nSubNodes(x),ncol=length(newcols)))
# 				names(newdat) <- newcols
# 				x@subnode.data = cbind(x@subnode.data,newdat)
# 			}
# 			
# 			# Remove superfluous columns
# 			if(length(supercols)!=0)
# 			{
# 				for(scol in supercols)
# 					x = rmdata(x,scol)
# 			}
# 		}else{
# 			# add this data:
# 			if(length(newcols)!=0)
# 			{
# 				print("adding this way")
# 				sndata(x,datnames=newcols) <- snode.data
# 			}
# 		}
# 		
# 		return(x)
# })
# 
# 
