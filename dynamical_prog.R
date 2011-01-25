require(phylobase)

# Implementation of dynamical programming algorithm (Sankoff, 1975)
# From algorithm laid out in Computational Molecular Evolution (Yang)
# See section 3.4 Maximum Parsimony of the above book for example details
#
# @author Conrad Stack
# @dedication: for practice...

dyn.loop.par<-function(phytree=junk,states.unique,states.freq,cost.mat,tip.states)
{
	phytree = as(phytree,'phylo4')
	# plot initial tree
	plot(phytree,show.node.label=T,rot=-90,tip.order=rev(seq(nTips(phytree))))
	phytree = reorder(phytree,"postorder")
	
	###############################################
	# Setting up data structures
	# --------------------------
	# "Cost" of each tip is simply it's 'from' row in the cost
	# matrix.  Or, simply, the respective cost of transition from this
	# one (and unchanging) state.
	#
	tip.costs = matrix(NA,nrow=nTips(phytree),ncol=length(states.unique))
	tip.states.mapping = sapply(tip.states,function(i) head(which(i==states.unique),1))
	for(ii in seq(nTips(phytree)))
	{
		tip.costs[ii,] = cost.mat[tip.states.mapping[ii],]
	}
	#tip.states.mat = matrix(rep(seq(length(states.unique)),nTips(phytree)),nrow=nTips(phytree),byrow=T)
	tip.states.mat = matrix(rep(states.unique,nTips(phytree)),nrow=nTips(phytree),byrow=T)
	rownames(tip.costs) <- tipLabels(phytree)
	rownames(tip.states.mat) <- tipLabels(phytree)
	
	# setup matrices to hold internal node costs and states
	#
	node.states = matrix(NA,nrow=nNodes(phytree),ncol=length(states.unique))
	node.costs = matrix(NA,nrow=nNodes(phytree),ncol=length(states.unique))
	rownames(node.costs) <- nodeLabels(phytree)
	rownames(node.states) <- nodeLabels(phytree)
	
	# taking advantage of phylo4 class: tips are numbered 1-nTips
	states = rbind(tip.states.mat,node.states)
	costs = rbind(tip.costs,node.costs)
	
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
		daughters = edges(phytree)[which(edges(phytree)[,1]==node.i),2]
		tmp.scores = numeric(length(states.unique))
		tmp.states = character(length(states.unique))
		
		cat("running on: ",labels(phytree)[node.i])	
		cat(" with daughters: ", labels(phytree)[daughters],"\n")
		
		# for internal nodes
		if(nodeType(phytree)[node.i] != "root")
		{
			# state.x is a potential mother state
			for (state.x in seq(length(states.unique)))
			{
				# state.y is potential state of node.i
				tmp.scores.y = numeric(length(states.unique))
				for(state.y in seq(length(states.unique)))
				{
					tmp.scores.y[state.y] = cost.mat[state.x,state.y]
					for(daught in daughters)
						tmp.scores.y[state.y] = tmp.scores.y[state.y] + costs[daught,state.y]
				}
	
				minimums = which(tmp.scores.y == min(tmp.scores.y))
				tmp.scores[state.x] = tmp.scores.y[which.min(tmp.scores.y)]
				tmp.states[state.x] = paste(states.unique[which(tmp.scores.y == min(tmp.scores.y))],collapse="/")
			}
		} else {
			
			# for root node:
			tmp.scores.y = numeric(length(states.unique))
			for(state.y in seq(length(states.unique)))
				for(daught in daughters)
					tmp.scores.y[state.y] = tmp.scores.y[state.y] + costs[daught,state.y]
			
			tmp.scores=tmp.scores.y
			tmp.states = states.unique
		}
		costs[node.i,] = tmp.scores
		states[node.i,] = tmp.states
	}
	
	return(list(states,costs))

}


### TEST dyn.loop

ntips = 6
junk = as(rcoal(ntips),'phylo4')
plot(junk,rot=-90)
nodeLabels(junk) <- paste("node",seq(from=nTips(junk)+1,to=nTips(junk)+nNodes(junk)),sep="")


# change to post-order
junk = reorder(junk,"postorder")
sunique=c("T","C","A","G")
sfreq = rep(1/4,4)
cmat = matrix(c(0,1,1.5,1.5,
			1,0,1.5,1.5,
			1.5,1.5,0,1,
			1.5,1.5,1,0),
			nrow=4,ncol=4,byrow=T)
rownames(cmat) <- paste("from",sunique,sep="-")
colnames(cmat) <- paste("to",sunique,sep="-")


# Choose random states
#tip.states = sample(states.unique,ntips,replace=T)
tstates = c("C","C","A","G","A","A") # from Yang, p. 997
opar <- par() 
par(las=2) # make label text perpendicular to axis
tipLabels(junk) <- paste(tipLabels(junk),tstates,sep=": ")
junk = read.tree("maxpars.tree")

dyn.loop.par(junk,sunique,sfreq,cmat,tstates)



# TODO:	-Change cost matrix to transition log-probability matrix
#		 see pagel 2004
#		-Each branch gets it's own cost matrix (weighted by the branch len)
#		-This is used to calculate the log-probability at each node.
#		-dyn.prog. algorithm then used to sum up most likely states
#		-MLE algorithm used to search state space?
#		-RE-READ pagel,pupko and try this approach
#

