require(phylobase)

rd = "(Taxon1:0.1, (Taxon2:0.3, Taxon3:0.4):0.5);"
new.tree = read.tree(text=rd)
plot(new.tree)

rd.simmap = "((Taxon1:{A,0.1; C,0.1}, Taxon2:{T,0.1; C,0.1}):{C,0.5}, Taxon3:{C,0.4} );"
sub.br = ":\\{([a-zA-z0-9]{1,5}),(\\d{0,10}\\.\\d{0,10});"
#gsub(rd.nowhite, sub.br ,"-(\\1)-(\\2)-" ,rd.nowhite)


#------------------------
# SIMMAP Methods		|
#------------------------

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
		treesblock = read.nexus.block(txt=rawtext,block="trees")
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


testtree = read.simmap(text="(((((32:{0,0.058907691963},(31:{0,0.022424212177:1,0.016767967599},(34:{1,0.036771632295},(30:{1,0.033928324112},33:{0,0.029725747655:1,0.004202576457}):{1,0.002843308183}):{1,0.002420547480}):{1,0.000275233450:0,0.019440278711}):{0,0.008552585573},29:{0,0.067460277504}):{0,0.040135207374},8:{0,0.091158418718:1,0.000880105051:0,0.015556961136}):{0,0.019704981181},(10:{0,0.112116752589},(2:{0,0.000772585273},3:{0,0.000772585273}):{0,0.014969340499:1,0.029030780722:0,0.067344046095}):{0,0.015183713480}):{0,0.023742583944},((1:{1,0.066050466671},4:{1,0.066050466671}):{1,0.039823924285},(9:{1,0.086482781129},(((5:{1,0.017325369522},6:{1,0.017325369522}):{1,0.006182126280},7:{1,0.023507495803}):{1,0.034790702136},(((26:{1,0.019800202825},14:{1,0.019800202825}):{1,0.016420404516},22:{1,0.036220607341}):{1,0.009769629593},(((23:{1,0.014983924906},(24:{1,0.004248281035},(18:{1,0.001668719645},21:{1,0.001668719645}):{1,0.002579561389}):{1,0.010735643875}):{1,0.009021350133},12:{1,0.024005275044}):{1,0.015472368742},(((((20:{1,0.019149815099},(15:{1,0.015578674920},11:{1,0.015578674920}):{1,0.003571140182}):{1,0.006407372276},25:{1,0.025557187375}):{1,0.001278231913},(13:{1,0.017974442632},16:{1,0.017974442632}):{1,0.008860976659}):{1,0.004074632240},(17:{1,0.027243985323},(28:{1,0.022983956116},19:{1,0.022983956116}):{1,0.004260029207}):{1,0.003666066205}):{1,0.005572724711},27:{1,0.036482776240}):{1,0.002994867546}):{1,0.006512593148}):{1,0.012307961005}):{1,0.028184583180}):{1,0.019391609822}):{1,0.019245403638:0,0.025923255434});")




#
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

	treesblock = read.nexus.block(txt=rawtext,block="trees")
	
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
			tipLabels(outtrees[[treeind]]) <- taxname[as.integer(tipLabels(outtrees[[treeind]]))]
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


# TODO: make this generic
# collapse singleton nodes
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



#-----------------------
# Process Nexus files  |
#-----------------------

# Method to read the first comment in a line in the format '[&...]'
# This is for reading tree weights chiefly
#
# Example:
# get.nexus.comments("example.txt")->lala
#
get.nexus.comments<-function(finput,text=NULL)
{
	if(!is.null(text)){
		# TODO: check text for newlines and split on them if they exist.
		rawtext=text
	} else {
		if(!file.exists(finput))
			stop("Assuming finput is a file and could not find it")
		
		rawtext = scan(finput,what=character(0),strip.white=T,sep="\n")		
	}
	
	# TODO: return named pair {treename, comment}
	comments = character(0)
	for(ii in seq(length(rawtext)))
	{
		junk =  gsub("^.*?\\[(.*?)\\].*$","\\1",rawtext[ii])
		if(length( grep("^&(.*)$",junk) ) != 0)
			comments = append(comments,junk)
	}
	
	return(comments)	
}


# Internal function - split tokens by a certain character
.split.tokens <- function(txt,char)
{
	nlines = length(txt)
	newtokens = character(nlines)
	charlines = unname(sapply(txt,function(i) length(grep(char,i))))
	curline=1
	for(ll in seq(nlines))
	{
		if(charlines[ll] == 0){
			newtokens[curline] = txt[ll]
			curline = curline + 1
		}else{
			tmp = strsplit(txt[ll],char)[[1]]
			for(tmpline in tmp){
				#if(tmpline!=""){
					newtokens[curline] = paste(tmpline,char,sep="")
					curline = curline + 1
				#}
			}
		}
	}
	return (newtokens)
}


# Internal function - get all content within a nexus block
.get.nexus.block.inds <- function(filename,blockname,text=NULL)
{
	# choose character vector
	if(!is.null(text))
	{
		filetext = text
	} else {
		filetext = scan(filename,what="character",sep="\n",strip.white=T)
	}
	filetext = tolower(filetext)
	
	start.ind = agrep(paste("begin",tolower(blockname)),filetext,ignore.case=T)
	if(length(start.ind) == 0)
		return (integer(0))
	
	end.ind = grep("end;",filetext,ignore.case=T)
	end.ind = end.ind[head(which(end.ind > start.ind),1)]

	return( c((start.ind),(end.ind)) )
}



# Internal function to parse / check assumptions block
.process.assumptions <- function(obj,block.txt)
{
	
	if(!is(obj[[1]],"brownie"))
		stop("Processed object needs to be of class brownie")
	
	
	for(aline in block.txt)
	{
		tokens = strsplit(aline,"\\s")[[1]]
		next.is.name=F
		next.is.taxa=F
		taxinds=numeric(0)
		taxname=""
		
		for(kk in seq(length(tokens)))
		{
			if(tolower(tokens[kk]) == "taxset")
			{
				next.is.name = T
				next
			} else {
				
				# reset and start recording names
				if(next.is.name){
					next.is.name=F
					taxname = sub("=","",tokens[kk])
					next.is.taxa=T
					next
				}
				
				if(next.is.taxa){
					if(length(grep(";",tokens[kk]))==1)
					{
						taxinds = append(taxinds,sub(";","",tokens[kk]))
						break
					} else {
						taxinds = append(taxinds,tokens[kk])
					}
				}
			}
		}
		
		# add this taxaset to the object, if there is one:
		if(taxname != "")
		{
			nm = paste("TAXSET",taxname,sep="_")
			taxaI = data.frame(all=rep(0,nTips(tree)))
			names(taxaI) <- nm
			
			for(tind in seq(length(obj)))
			{
				# convert if needed
				if(!inherits(obj[[tind]],"phylo4d"))
					obj[[tind]] = phylo4d(obj[[tind]])
				
				taxaI[,1] = sapply(tipLabels(obj[[tind]]),function(i) ifelse(i %in% taxinds,1,0),simplify=T)
				#names(taxaI) <- nm
				obj[[tind]] = addData(obj[[tind]],tip.data=taxaI)
				
			}
		}
	}
	
	return (obj)
}




# internal function to parse / check brownie block
.process.brownie <- function(obj,block.txt)
{
	return(obj)
}


# Internal:
# Get text of a read nexus block
read.nexus.block<-function(finput,txt=NULL,block)
{	
	if(!is.null(txt)){
		# Using the text argument is not recommended
		rawtext=txt
	} else {
		
		if(!file.exists(finput))
			stop("Assuming finput is a file and could not find it")
		
		rawtext = scan(finput,what=character(0),strip.white=T,sep="\n")		
	}
	
	inds = .get.nexus.block.inds(blockname=block,text=rawtext)
	if(length(inds)==0)
	{
		warning(paste("This file has",block, "no block"))
		return (character(0))
	}
	
	# TODO: split up newlines if they exist
	return (rawtext[(inds[1]+1):(inds[2]-1)])
}






#-----------------------
# Read / Write brownie |
#-----------------------


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
		
		# should return a list
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
	brau.new = .process.assumptions(brau.new,assumptions.part)  # TODO: overload constructor to do this processing.
	
	return(brau.new)
}


