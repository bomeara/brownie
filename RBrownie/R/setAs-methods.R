#------------------------
# Coersing functions	|
#------------------------

### to/from brownie
#
setAs('brownie','phylo',function(from,to) {
	
	# brownie will be automatically coersed to phylo4d,
	# then conversion functions are handled by phylobase 
	# methods
	#
	phyobj = as(from,'phylo4')
	apeobj = NULL
	
	# TODO: Should this collapse singles?  Or create zero-length branches?
	if(hasSingle(phyobj)){
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
	#res@annote <- annote
	
	return(res)
})


## from/to phylo4d_ext:
#
# up convert:


# down convert:
# setAs('phylo4d_ext','phylo4d',function(from,to) {
# 	
# 	#TODO: Need to convert subnodes back into singletons
# 
# })
# 
# 
# setAs('phylo4d_ext','phylo4',function(from,to) {
# 	
# 	#TODO: Need to convert subnodes back into singletons
# 
# })
# 
# 
# setAs('phylo4d_ext','phylo',function(from,to) {
# 	
# 	# TODO: Need to convert subnodes back into singlestons, 
# 	#		and then collapse the singletons.
# 	
# })

