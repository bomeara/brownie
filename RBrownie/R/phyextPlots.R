#----------------------------------------
#  Phylo4d extension plots
#----------------------------------------

# NOTE: datapart can be a column index(?) or character string
#
# @param datapart can be any integer specifying the column index of the data or a character string 
#
# TODO: use par(bg = "color") to discover background color of plotting device
#
phyextPlot <- function(x,states,states.col,
						states.na="none", 
						datapart=1,
						plot.subnodes=T,
						plot.points=T,
						line.widths,line.types, ... )
{

	# plot base phylogeny using phylobase functions:
	junk <- x
		
	gtree <- extractTree(junk)  # TODO: overload this so that it works on phyext class (keeps subnode info)
	posi = phyloXXYY(gtree)
	grid.newpage()
	pushViewport(viewport(layout=grid.layout(nrow=1, ncol=1), name="base"))
	pushViewport(viewport(layout.pos.col=1, name="plot1"))
		
	treePlot(gtree, newpage=FALSE,...)

	if(hasTipData(junk) || hasNodeData(junk))
	{
		if(missing(states))
		{
			tmp = as.character(unique(tdata(junk,'all')[,datapart,drop=F])[,1])
			tmp = tmp[!is.na(tmp)]
			states = c(states.na,tmp)
			states.col = c(1,seq(from=2,length.out=length(states)-1))
		}
		
		if(missing(line.widths)){
			line.widths = rep(1,length(states))
		} else {
			if(length(states) != length(line.widths))
				stop("line.widths need to be the same length as states")
		}
		
		if(missing(line.types)){
			line.types = rep(1,length(states))
		} else {
			if(length(states) != length(line.types))
				stop("line.types need to be the same length as states")
		}
				
		seekViewport("tree")
		
		eord = edges(gtree)[posi$eorder,] # this is the order used
		treedata = tdata(junk,"all")[eord[,2],datapart,drop=T]  
		datamap = sapply(treedata,function(i) which(states == i),simplify=T)
		if(is.list(datamap))
			datamap = unlist(lapply(datamap, function(i) ifelse(length(i)==0,1,i[1])))
		
		# replot edges:
		grid.segments(posi$segs$h0x,posi$segs$h0y, posi$segs$h1x,posi$segs$h1y,gp=gpar(col=states.col[datamap],lwd=line.widths[datamap],lty=line.types[datamap]))
		grid.segments(posi$segs$v0x,posi$segs$v0y, posi$segs$v1x,posi$segs$v1y,gp=gpar(col=states.col[datamap],lwd=line.widths[datamap],lty=line.types[datamap]))
		if(plot.points) grid.points(posi$xx,posi$yy,pch=20,gp=gpar(col=states.col[datamap],cex=0.5))
		
		
		# plot sub nodes:
		if(hasSubNodes(junk) && plot.subnodes)
		{
			esub = edges(junk)[getSubNodeEdgeInds(junk),]
			posi.inds = apply(esub,1,function(i) which(i[1] == eord[,1] & i[2] == eord[,2]))
			subdata = getSubNodeData(junk,datapart)
			submapping = sapply(subdata[,1],function(i) which(states == i))
			subposi = getSubNodePosish(junk)
			
			get.x.offset <- function(xxyy,inds)
			{
				apply(cbind(xxyy$segs$h0x[inds],xxyy$segs$h1x[inds]),1,diff)
			}
			
			# reorder (so that lines don't completely cover each other):
			neword = order(rowMeans(subposi),decreasing=T)
			subposi = subposi[neword,]
			submapping = submapping[neword]
			posi.inds = posi.inds[neword]
			
			# subbranch positions:
			subposi.x0 = posi$segs$h0x[posi.inds]
			subposi.y0 = posi$segs$h0y[posi.inds]
			subposi.y1 = posi$segs$h1y[posi.inds]
			subposi.x1 = posi$segs$h0x[posi.inds] + (get.x.offset(posi,posi.inds) * rowMeans(subposi))
			subposi.vx0 = posi$segs$v0x[posi.inds]
			subposi.vy0 = posi$segs$v0y[posi.inds]
			subposi.vy1 = posi$segs$v1y[posi.inds]
			subposi.vx1 = posi$segs$v1x[posi.inds]
			
			# plot subnodes:
			grid.segments(subposi.x0,subposi.y0,subposi.x1,subposi.y1,gp=gpar(col=states.col[submapping],lwd=line.widths[submapping],lty=line.types[submapping]))
			grid.segments(subposi.vx0,subposi.vy0,subposi.vx1,subposi.vy1,gp=gpar(col=states.col[submapping],lwd=line.widths[submapping],lty=line.types[submapping]))
			if(plot.points) grid.points(subposi.x1,subposi.y1,pch=20,gp=gpar(col=states.col[submapping],cex=0.5))
			
		}
		
		upViewport(2)
	}
}


setGeneric('plot')
setMethod('plot', signature(x='phylo4d_ext', y='missing'), function(x, y, ...) {
    phyextPlot(x, ...)
})


plot.taxaset <- function(x,taxind,taxcol="red",taxlwd=1,excol="grey",exlwd=1,blankit=F,...)
{
	if(hasTaxasets(x))
	{
		index = taxind.to.dataind(x=x,taxind=taxind)
		
		# Rename internal nodes:
		#cat("Renaming internal nodes (might take a while if there are a lot of taxa):\n")
		###########################
		taxmrca = integer(0)
		taxnames = taxa.charvect(x,taxind)
		top = unname(MRCA(x,taxnames))
		if(!areTaxaMono(x,taxind)){
			taxout = setdiff(names(descendants(x,top)),taxnames)
			taxmrca = unname(ancestor(x,MRCA(x,taxout)))
		}
		
		taxinds = taxaname.to.taxind(x,taxnames)
		lambdaS = function(nodetwo,tree,nodeone) shortestPath(tree,nodeone,nodetwo)
		taxmrca = c(taxmrca,unique(unlist(sapply(taxinds,lambdaS,x,top))))
		taxmrca = c(taxmrca,taxinds)

		if(top %in% taxmrca)
			taxmrca = taxmrca[-which(taxmrca==top)]
		tdata(x)[,index][is.na(tdata(x)[,index])] = 0
		tdata(x)[,index][taxmrca] = 1 
		#####################
		
		if(!blankit){
			phyextPlot(x,states=c(0,1),states.col=c(excol,taxcol),datapart=index,plot.subnodes=F,line.widths=c(exlwd,taxlwd),plot.points=F,...)
		} else {
			nodecutcol=par("bg")
			if(nodecutcol=="transparent") nodecutcol="white"
			tdata(x)[,index][top] = -1
			phyextPlot(x,states=c(0,1,-1),states.col=c(excol,taxcol,nodecutcol),datapart=index,plot.subnodes=F,line.widths=c(exlwd,taxlwd,1),plot.points=F,...)
		}
		grid.text(paste("Taxaset:",sprintf('%s',taxaset.names(x)[taxind])),
					x=unit(2, "mm"), y=unit(1, "npc") - unit(2, "mm"),
           			just=c("left", "top"))
	} else {
		warning("Tree x does not contain taxasets (is it a brownie object?).")
	}
}

#
#
plot.censored <- function(x,taxind,
							taxcol="red",
							taxlwd=1,
							excol="grey",
							exlwd=1,
							rmlty=3,
							...)
{
	if(hasTaxasets(x))
	{
		showcut=F
		index = taxind.to.dataind(x=x,taxind=taxind)
		
		# Rename internal nodes:
		#cat("Renaming internal nodes (might take a while if there are a lot of taxa):\n")
		###########################
		taxmrca = integer(0)
		taxnames = taxa.charvect(x,taxind)
		top = unname(MRCA(x,taxnames))
		if(areTaxaPara(x,taxind)){
			
			# TODO: make sure subtop/other top covers the range of missing nodes.
			taxout = setdiff(names(descendants(x,top)),taxnames)
			taxout.mrca = MRCA(x,taxout)
			taxmrca = ancestor(x,taxout.mrca)
			subtop = taxmrca
			other.top = setdiff(children(x,taxmrca),taxout.mrca)  # large monophyly
			showcut=T
			
		}
				
		taxinds = taxaname.to.taxind(x,taxnames)
		lambdaS = function(nodetwo,tree,nodeone) shortestPath(tree,nodeone,nodetwo)
		taxmrca = c(taxmrca,unique(unlist(sapply(taxinds,lambdaS,x,top))))
		taxmrca = c(taxmrca,taxinds)

		if(top %in% taxmrca)
			taxmrca = taxmrca[-which(taxmrca==top)]
		
		tdata(x)[,index][is.na(tdata(x)[,index])] = 0
		tdata(x)[,index][taxmrca] = 1 
		#####################
		
		nodecutcol=par("bg")
		if(nodecutcol=="transparent") nodecutcol="white"
		
		if(showcut)
		{
			tdata(x)[,index][subtop] = -1
			tdata(x)[,index][other.top] = -1
			tdata(x)[,index][taxout.mrca] = -1
		}

		phyextPlot(x,states=c(0,1,-1),
					states.col=c(excol,taxcol,nodecutcol),
					datapart=index,
					plot.subnodes=F,
					line.widths=c(exlwd,taxlwd,1),
					line.types=c(1,1,rmlty),
					plot.points=F,
					edge.color=excol,...)
		
		grid.text(paste("Taxaset:",sprintf('%s',taxaset.names(x)[taxind])),
					x=unit(2, "mm"), y=unit(1, "npc") - unit(2, "mm"),
           			just=c("left", "top"))
	} else {
		warning("Tree x does not contain taxasets (is it a brownie object?).")
	}
}

