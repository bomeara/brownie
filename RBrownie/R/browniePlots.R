#----------------------------------------
#  brownie plots
#----------------------------------------


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

