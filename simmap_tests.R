require(phylobase)

rd = "(Taxon1:0.1, (Taxon2:0.3, Taxon3:0.4):0.5);"
new.tree = read.tree(text=rd)
plot(new.tree)

rd.simmap = "((Taxon1:{A,0.1; C,0.1}, Taxon2:{T,0.1; C,0.1}):{C,0.5}, Taxon3:{C,0.4} );"

sub.br = ":\\{([a-zA-z0-9]{1,5}),(\\d{0,10}\\.\\d{0,10});"
#gsub(rd.nowhite, sub.br ,"-(\\1)-(\\2)-" ,rd.nowhite)

# get node names from newick string
get.nodenames<-function(newick.txt)
{
	nnames = character(0)
	ttmp = newick.txt
	nname.pat="(\\(|,|\\))([a-zA-Z0-9]{1,}):"
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


# read modified newick file
# citation: 
# Bollback J. P. (2006) SIMMAP: Stochastic character mapping of 
# discrete traits on phylogenies. BMC Bioinformatics. 7:88
# 
# Assume that branch lengths are present (otherwise, why use SIMMAP?)
#
read.simmap <- function(file="",text=NULL, version=1.1, ...)
{
	if(is.null(text))
		stop("Need to have text for now")

	# clear whitespace	
	text = gsub("\\s","",text)
	
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
		splitchar = ifelse(version==1.0,";",":")
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
	#
	##
	
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
	nnames = head(get.nodenames(newick.str),-1)
	scount = 1
	nTips = length(tr$tip.label)
	nInts = length(tr$node.label)

	# loop through split branches
	for(kk in unique(edge.ind))
	{
		# this information shouldn't change:
		edgeind = which(edge.ind == kk)
		trind = which(nnames[kk] == all.labels)
		esplice = which(tr$edge[,2] == trind) # except this, this will change
		anc = tr$edge[esplice,1]
		dec = tr$edge[esplice,2]
	
		newsingles = length(edgeind)-1
		if(newsingles==0)
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
		
			# update info:
			esplice = length(tr$edge.length) # should be the last one
			dec = nodeid
			scount = scount + 1
		}
	}
	
	write.nexus(tr,file="written.tree")
	return(as(tr,"phylo4"))
}

testtree = read.simmap(text="(((((32:{0,0.058907691963},(31:{0,0.022424212177:1,0.016767967599},(34:{1,0.036771632295},(30:{1,0.033928324112},33:{0,0.029725747655:1,0.004202576457}):{1,0.002843308183}):{1,0.002420547480}):{1,0.000275233450:0,0.019440278711}):{0,0.008552585573},29:{0,0.067460277504}):{0,0.040135207374},8:{0,0.091158418718:1,0.000880105051:0,0.015556961136}):{0,0.019704981181},(10:{0,0.112116752589},(2:{0,0.000772585273},3:{0,0.000772585273}):{0,0.014969340499:1,0.029030780722:0,0.067344046095}):{0,0.015183713480}):{0,0.023742583944},((1:{1,0.066050466671},4:{1,0.066050466671}):{1,0.039823924285},(9:{1,0.086482781129},(((5:{1,0.017325369522},6:{1,0.017325369522}):{1,0.006182126280},7:{1,0.023507495803}):{1,0.034790702136},(((26:{1,0.019800202825},14:{1,0.019800202825}):{1,0.016420404516},22:{1,0.036220607341}):{1,0.009769629593},(((23:{1,0.014983924906},(24:{1,0.004248281035},(18:{1,0.001668719645},21:{1,0.001668719645}):{1,0.002579561389}):{1,0.010735643875}):{1,0.009021350133},12:{1,0.024005275044}):{1,0.015472368742},(((((20:{1,0.019149815099},(15:{1,0.015578674920},11:{1,0.015578674920}):{1,0.003571140182}):{1,0.006407372276},25:{1,0.025557187375}):{1,0.001278231913},(13:{1,0.017974442632},16:{1,0.017974442632}):{1,0.008860976659}):{1,0.004074632240},(17:{1,0.027243985323},(28:{1,0.022983956116},19:{1,0.022983956116}):{1,0.004260029207}):{1,0.003666066205}):{1,0.005572724711},27:{1,0.036482776240}):{1,0.002994867546}):{1,0.006512593148}):{1,0.012307961005}):{1,0.028184583180}):{1,0.019391609822}):{1,0.019245403638:0,0.025923255434});")



# expand singleton nodes into bifurcating nodes with one junk node
# this function is mainly for plotting and saving files
# TODO: convert all explicit calls to @'slot' to their abstract counterpart (e.g. tree@edge.length goes to edgeLength(tree)
#
expand.singles <- function(tree)
{
	# note: tips should be indexed 1...N, where N is the number of tips
	if(!is(tree,"phylo4"))
		stop("tree needs to be of class phylo4")
	
	tmptable=table(tree@edge[,1])
	snodes = as.integer(names(tmptable)[which(tmptable==1)])
	snodes = snodes[snodes!=0] # 0 is the root node (check on this...)
	
	nold = nrow(tree@edge)
	nnew = length(snodes)
	count = 1
	for(ii in snodes)
	{
		tree@label = append(sprintf("Junk%0.7d",count),tree@label)
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
	names(tree@label) <- as.character(seq(1,nnodes))
	names(tree@edge.length) <- apply(tree@edge,1,paste,collapse="-")
	
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




#-----------------------
# Read from nexus file |
#-----------------------

# Method to read the first comment in a line in the format '[&...]'
# This is for reading tree weights chiefly
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
#get.nexus.comments("example.txt")->lala  # test call


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


# internal function to parse / check assumptions block
.process.assumptions <- function(obj,block.txt)
{
	if(!is(obj,"brownie"))
		stop("Processed object needs to be of class brownie")
	
	return (obj)
}

# internal function to parse / check brownie block
.process.brownie <- function(obj,block.txt)
{
	return(obj)
}


# read nexus assumptions block
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


read.brownie<-function(fname)
{	
	if(!file.exists(fname))
		stop(paste("File",fname,"cannot be found in",getwd()))
	
	filetxt = scan(fname,what=character(0),strip.white=T,sep="\n")
	brownie.part = read.nexus.block(txt=filetxt,block="BROWNIE")
	assumptions.part = read.nexus.block(txt=filetxt,block="ASSUMPTIONS")
	
	phy.part = readNexus(fname)  # convert nexus strings to simmap tree
	
	brau.new = new("brownie",phy.part,commands=brownie.part)
	brau.new = .process.assumptions(assumptions.part)  # TODO: overload constructor to do this processing.
	
	return(brau.new)
}


