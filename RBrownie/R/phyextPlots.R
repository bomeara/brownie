#----------------------------------------
#  Phylo4d extension plots
#----------------------------------------

# 
# @param usestate specifies which column index for the @data slot should be accessed and plotted.  
#		 It can be either an integer (between 1...ncol) or a character string (which is a column 
#		 name.)
#
phyextPlot <- function(x,states,states.col,
						states.na="none", 
						usestate=1,
						plot.subnodes=T,
						plot.points=T,
						line.widths,line.types, ... )
{

	# TODO: sanity check to make sure usestate is a proper data column index
	
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
			tmp = as.character(unique(tdata(junk,'all')[,usestate,drop=F])[,1])
			tmp = tmp[!is.na(tmp)]
			if(hasSubNodes(junk) && plot.subnodes)
			{
				tmp = c(tmp,as.character(unique(sndata(junk)[,usestate,drop=F])[,1]))
				tmp = tmp[!is.na(tmp)]
			}
			tmp = unique(tmp)
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
		treedata = tdata(junk,"all")[eord[,2],usestate,drop=T]  
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
			esub = matrix(edges(junk)[getSubNodeEdgeInds(junk),],ncol=2)
			posi.inds = apply(esub,1,function(i) which(i[1] == eord[,1] & i[2] == eord[,2]))
			subdata = getSubNodeData(junk,usestate)
			submapping = sapply(subdata[,1],function(i) which(states == i))
			subposi = getSubNodePosish(junk)
			
			get.x.offset <- function(xxyy,inds)
			{
				apply(cbind(xxyy$segs$h0x[inds],xxyy$segs$h1x[inds]),1,diff)
			}
			
			# reorder (so that lines don't completely cover each other):
			neword = order(rowMeans(subposi),decreasing=T)
			subposi = matrix(subposi[neword,],ncol=2)
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

# The line below is unnecessary since phylobase sets 'plot' as a generic function
# (through graphics I think)
#setGeneric('plot')
setMethod('plot', signature(x='phylo4d_ext', y='missing'), function(x, y, ...) {
    phyextPlot(x, ...)
})

