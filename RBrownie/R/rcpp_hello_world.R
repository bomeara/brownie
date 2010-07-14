
rcpp_hello_world <- function(){
	.Call( "rcpp_hello_world", PACKAGE = "RBrownie" )
}

read.brownie <- function(filename)
{
	# check argument for character string
	if( !is.character(filename) && length(filename) == 1)
		stop(paste("first argument must be a character string specifying one file only"))
	
	# check if filename exists
	if( !file.exists(filename) )
		stop( paste("filename,",filename,", does not exist") )
	
	# convert to unix-style path
	filename = gsub("\\\\","/",filename)
	
	# add single quotes (brownie requirement for full paths
	filename = paste("'",filename,"'",sep="")
	
	.Call("readBrownie",filename,PACKAGE="RBrownie")
}
