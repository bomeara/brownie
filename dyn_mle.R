# Dynamic programming: infer maximum likelihood
#
# See Computational Molecular Evolution (Yang)
# 	Section 4 - Maximum Likelihood Methods
require(phylobase)

#
# ripped from geiger (thanks luke)
logspace_add<-function(logx, logy) {
	if(logx==-Inf) return(logy) else max(logx, logy) + log1p(exp (-abs (logx - logy)));
}

logspace_sum<-function(logx) {
      r<-logx[1]
      if(length(logx)>1)
      	for(i in 2:length(logx))
      		r<-logspace_add(r, logx[i])
      r	
}


dyn.loop.mle<-function(phytree, cost.mat, states.unique, states.freq, tip.states=NULL,conditionals=T)
{
	
	phytree = as(phytree,'phylo4')
	# plot initial tree
	#plot(phytree,show.node.label=T,rot=-90,tip.order=rev(seq(nTips(phytree))))
	phytree = reorder(phytree,"postorder")


	###############################################
	# Setting up data structures
	# --------------------------
	# "Cost" of each tip is simply it's 'from' row in the cost
	# matrix.  Or, simply, the respective cost of transition from this
	# one (and unchanging) state.
	#
	
	tip.costs = matrix(NA,nrow=nTips(phytree),ncol=length(states.unique))
	tip.states.mapping = sapply(tip.states,function(i) head(which(i==states.unique),1),simplify=T)
		
	cat("mapping states: ",tip.states.mapping,"\n")
	for(ii in seq(nTips(phytree)))
	{
		tip.costs[ii,] = 0
		tip.costs[ii,tip.states.mapping[ii]] = 1
	}
	tip.states.mat = matrix(rep(states.unique,nTips(phytree)),nrow=nTips(phytree),byrow=T)
	rownames(tip.costs) <- nodeId(phytree,"tip")
	colnames(tip.costs) <- states.unique
	


	# setup matrices to hold internal node costs and states
	#
	#node.states = matrix(NA,nrow=nNodes(phytree),ncol=length(states.unique))
	node.costs = matrix(NA,nrow=nNodes(phytree),ncol=length(states.unique))
	rownames(node.costs) <- nodeId(phytree,"internal")
	#rownames(node.states) <- nodeLabels(phytree)
	
	# taking advantage of phylo4 class: tips are numbered 1-nTips
	#states = rbind(tip.states.mat,node.states)
	costs = rbind(tip.costs,node.costs)
	colnames(costs) <- states.unique

	####################################
	# Dynamical programming algorithm
	# -------------------------------
	#
	#
	# for each edge:
	for(jj in seq(nEdges(phytree)))
	{
		node.i = edges(phytree)[jj,2]

		# already calculated tip scores when 
		# data structures were set up
		if(nodeType(phytree)[node.i] == "tip")
			next
				
		# temporary data structures
		edge.inds = which(edges(phytree)[,1] == node.i)
		daughters = edges(phytree)[edge.inds,2]
		tmp.scores = numeric(length(states.unique))
		#tmp.states = character(length(states.unique))
		
		cat("running on: ",labels(phytree)[node.i])	
		cat(" with daughters: ", labels(phytree)[daughters],"\n")
		
		tmp.costs = numeric(length(states.unique))			
		tmp.costmat = array(NA,c(length(edge.inds),dim(cost.mat)))
		
		# generate temporary transition probability matrix for these branches (log-like)
		#
		for(kk in seq(length(edge.inds)))
		{
			elen = edgeLength(phytree)[edge.inds[kk]]
			tmp.costmat[kk,,] = matexpo( cost.mat * elen) 
		}
		
		# state.x is a potential mother state
		# state.y is a potential child state
		for (state.x in seq(length(states.unique)))
		{
			tmp.score.x = 1
			for(kk in seq(length(edge.inds))){
				tmp.scores = 0
				for(state.y in seq(length(states.unique)))
					tmp.scores = tmp.scores + (tmp.costmat[kk,state.y,state.x] * costs[daughters[kk],state.y])
				
				tmp.score.x = tmp.score.x * tmp.scores
			}
			tmp.costs[state.x] = tmp.score.x
		}

		# store this nodes costs:
		#costind = nTips(phytree) + which(node.i == nodeId(phytree))
		costs[node.i,] = tmp.costs
	}
	

	lik = sum( costs[rootNode(phytree),] * states.freq )
	# reconstruct vector:
	#outframe = data.frame(unname(cbind(costs,states)),stringsAsFactors=F)
	#colnames(outframe) <- c(paste("loglike",states.unique,sep="."),paste("state",states.unique,sep="."))
	#return( addData(phytree,all.data = outframe) )
	if(conditionals)
		return(costs)
	
	return(lik)
}


##### MLE TEST

# Tree:
treestr = "(((1:0.2,2:0.2):0.1,3:0.2):0.1,(4:0.2,5:0.2):0.1);"
tree = read.tree(text=treestr)
treed = phylo4d(tree,tip.data=data.frame(tipdata=c("T","C","A","C","C")))
nodeLabels(treed)[1] <- 0
nodeLabels(treed)[2:nNodes(treed)] <- nTips(treed) + (2:nNodes(treed)) - 1


# DNA Evoluation model:
# K80
states.freq = rep(1/4,4)
states.unique=c("T","C","A","G")
k = 2
rateq = matrix(c(0,k,1,1,
			k,0,1,1,
			1,1,0,k,
			1,1,k,0),nrow=4,byrow=T)

diag(rateq) = -apply(rateq,1,sum)
rateq = rateq %*% diag(states.freq)
p01 = matexpo(rateq*0.1)
p02 = matexpo(rateq*0.2)
p01
p02

tipLabels(treed) <- paste(tipLabels(treed),as.character(treed@data$"tipdata"),sep=": ")

plot(as(treed,'phylo4'),show.node.label=T)

dyn.loop.mle(treed,rateq,states.unique,states.freq,as.character(treed@data$"tipdata"))


#
#
#
#
#
#
#
###
##### Test: example from pupko 2000
###
#
#
#

states.unique = c("A","V")
states.freq = c(0.6,0.4)
cmat = matrix(c(0.7,0.3,0.45,0.55),ncol=2,byrow=T)

tip.states = c("V","A","V","A")
textstr = "(2:1,(3:1,(4:1,5:1):1):1);"
tree=read.tree(text=textstr)
#tree=root(tree,node=6,resolve.root=T)
phytree=as(tree,'phylo4')
nodeLabels(phytree) <- as.character(seq(nTips(phytree)+2,nNodes(phytree)+nTips(phytree)+1))
tipLabels(phytree) <- paste(tipLabels(phytree),tip.states)

plot(phytree,show.node.label=T)


dyn.loop.mle(phytree,cmat,states.unique,states.freq,tip.states)->junk

tree <- reorder(junk,"preorder")
ntypes = nodeType(phytree)
all.states = c(tip.states,character(nNodes(phytree)))
nstates = length(states.unique)

for(jj in seq(nrow(edges(tree))) )
{
	anc=tree@edge[jj,1]
	dec=tree@edge[jj,2]
	
	# don't worry about tips
	if(ntypes[dec]=="tip")
		next
	
	if(anc == 0)
	{
		# root
		llikes = as.numeric(tree@data[dec,seq(nstates)])
		chars = as.character(tree@data[dec,seq(nstates+1,2*nstates)])
		all.states[dec] = chars[which.max(llikes)]
	} else {
		all.states[dec] = as.character(tree@data[dec,seq(nstates+1,2*nstates)])[which(all.states[anc] == states.unique)]
	}	
}



