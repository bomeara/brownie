#------------------------
# Exporting functions	|
#------------------------

setAs('brownie','phylo',function(from,to) {
	
	# brownie will be automatically coersed to phylo4d,
	# then conversion functions are handled by phylobase 
	# methods
	#
	phyobj = as(from,'phylo4')
	apeobj = NULL
	
	# TODO: Should this collapse singles?  Or create zero-length branches?
	if(hasSingle(junk)){
		#apeobj = collapse.singles(apeobj)
		apeobj = as(expand.singles(phyobj),'phylo')
	} else {
		apeobj = as(phyobj,'phylo')
	}
	apeobj
})



## first arg is a phylo
setAs('phylo','brownie',function(from,to) {
	
	# this method will also strip labels which start with JUNK
	has.junk = any(as.logical(sapply(from$tip.label, function(i) length(grep("^JUNK.*",i)), simplify=T )))
	if(any(from$edge.length==0) || has.junk){
		res1 <- collapse.to.singles(x)
		res2 <- as(res1, "phylo4d")
		res <- new("brownie",res2,commands=character(0),datatypes=character(0))
	} else {
		res <- new("brownie",as(x,'phylo4d'),commands=character(0),datatypes=character(0))
	}
	
	#TODO?: make default annote arg NULL, and only assign if !is.null;
	# then update phylo4d methods accordingly (same thing with metadata?)
	res@annote <- annote
	
	return(res)
})



