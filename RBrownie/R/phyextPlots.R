#----------------------------------------
#  Phylo4d extension plots
#----------------------------------------

# NOTE: datapart can be a column index(?) or character string
#
phyextPlot <- function(x,states,states.col,states.na="none", datapart='simmap_state', ... )
{

	# plot base phylogeny using phylobase functions:
	junk <- x
	
	if(missing(states))
	{
		tmp = as.character(unique(tdata(junk,'all')[datapart])[,1])
		tmp = tmp[!is.na(tmp)]
		states = c(states.na,tmp)
		states.col = c(1,seq(from=2,length.out=length(states)-1))
	}
	
	
	gtree <- extractTree(junk)  # TODO: overload this so that it works on phyext class (keeps subnode info)
	posi = phyloXXYY(gtree)
	grid.newpage()
	pushViewport(viewport(layout=grid.layout(nrow=1, ncol=1), name="base"))
	pushViewport(viewport(layout.pos.col=1, name="plot1"))	
	treePlot(gtree, newpage=FALSE, ...)

	seekViewport("tree")
	
	eord = edges(gtree)[posi$eorder,] # this is the order used
	treedata = tdata(junk,"all")[eord[,2],datapart]  
	datamap = sapply(treedata,function(i) which(states == i),simplify=T)
	if(is.list(datamap))
		datamap = unlist(lapply(datamap, function(i) ifelse(length(i)==0,1,i[1])))
	
	# replot edges:
	grid.segments(posi$segs$h0x,posi$segs$h0y, posi$segs$h1x,posi$segs$h1y,gp=gpar(col=states.col[datamap]))
	grid.segments(posi$segs$v0x,posi$segs$v0y, posi$segs$v1x,posi$segs$v1y,gp=gpar(col=states.col[datamap]))
	grid.points(posi$xx,posi$yy,pch=20,gp=gpar(col=states.col[datamap],cex=0.5))
	
	
	# plot sub nodes:
	if(hasSubNodes(junk))
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
		grid.segments(subposi.x0,subposi.y0,subposi.x1,subposi.y1,gp=gpar(col=states.col[submapping]))
		grid.segments(subposi.vx0,subposi.vy0,subposi.vx1,subposi.vy1,gp=gpar(col=states.col[submapping]))
		grid.points(subposi.x1,subposi.y1,pch=20,gp=gpar(col=states.col[submapping],cex=0.5))
		
		
		upViewport(2)	
	}
		
}


setGeneric('plot')
setMethod('plot', signature(x='phylo4d_ext', y='missing'), function(x, y, ...) {
    phyextPlot(x, ...)
})

