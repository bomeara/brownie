#------------------------
# Coersing functions	|
#------------------------

# pulled from implicit(?) 'coerce<-' function
deep.phy.copy <- function (from, to, value,usedata=TRUE) 
{
	useslots = c("data", "metadata", "edge", "edge.length","label", "edge.label", "annote")	 # leaving 'order' out...
	if(!usedata) 
		useslots = useslots[-c(1,2)]
    for (what in useslots) 
    	slot(from, what) <- slot(value, what)
    from
}

# up convert:
# NOTE: not sure if it is necessary to define all of these...
# 		there might be a more 'S4' way to do this.
#
setAs('phylo','phylo4d_ext', function(from, to) {
	return(phyext(from))
})

setAs('phylo4','phylo4d_ext', function(from, to) {
	return(phyext(from))
})

setAs('phylo4d','phylo4d_ext', function(from, to) {
	return(phyext(from))
})

# 
# # down convert:
# # TODO: figure out if these functions are necessary:
# #
# setAs('phylo4d_ext','phylo4d',function(from,to) {
# 
# 	newphy = deep.phy.copy(new(to,order="unknown"),to,from,TRUE)
# 	newphy
# })
# 
# 
# setAs('phylo4d_ext','phylo4',function(from,to) {
# 	
# 	newphy = deep.phy.copy(new(to,order="unknown"),to,from,FALSE)
# 	newphy
# })
# 
# 
# 
# setAs('phylo4d_ext','phylo',function(from,to) {
# 	
# 	# TODO: Need to convert subnodes back into singlestons, 
# 	#		and then collapse the singletons.
# 	return(as(as(from,"phylo4d"),"phylo"))
# })
# 
# 




