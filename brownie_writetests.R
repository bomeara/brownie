require(RBrownie)
nexfiles = c("geospiza.nex","parrot.nex")
outfile = "write_tests.log"
cat("Starting write tests....\n\n",file=outfile)

for(ff in nexfiles)
{
	cat("\n",ff,"\n",file=outfile,append=T)
	test1 = readBrownie(ff)
	
	writeBrownie(test1,file=paste(ff,".tmp",sep=""))
	test1.remix = readBrownie(paste(ff,".tmp",sep=""))
	
	if(is.list(test1))
	{
		cat("--",sort(names(tdata(test1[[1]]))),"\n",file=outfile,append=T)
		cat("--",sort(names(tdata(test1.remix[[1]]))),"\n",file=outfile,append=T)
		cat("----------------------------------\n",file=outfile,append=T)
		cat("datatypes: ",datatypes(test1[[1]]),"\n",file=outfile,append=T)
		cat("datatypes: ",datatypes(test1.remix[[1]]),"\n",file=outfile,append=T)
		cat("----------------------------------\n",file=outfile,append=T)
		cat("commands: ",commands(test1[[1]]),"\n",file=outfile,append=T)
		cat("commands: ",commands(test1.remix[[1]]),"\n",file=outfile,append=T)
		cat("----------------------------------\n",file=outfile,append=T)
		cat("-tips ....",all(tipLabels(test1[[1]]) == tipLabels(test1.remix[[1]])),"\n")
		cat("----------------------------------\n",file=outfile,append=T)
		cat("----------------------------------\n",file=outfile,append=T)
		cat("weights: ",weight(test1),"\n",file=outfile,append=T)
		cat("weights: ",weight(test1.remix),"\n",file=outfile,append=T)	
	}else{
		cat("--",sort(names(tdata(test1))),"\n",file=outfile,append=T)
		cat("--",sort(names(tdata(test1.remix))),"\n",file=outfile,append=T)
		cat("----------------------------------\n",file=outfile,append=T)
		cat("datatypes: ",datatypes(test1),"\n",file=outfile,append=T)
		cat("datatypes: ",datatypes(test1.remix),"\n",file=outfile,append=T)	
		cat("----------------------------------\n",file=outfile,append=T)
		cat("commands: ",commands(test1),"\n",file=outfile,append=T)
		cat("commands: ",commands(test1.remix),"\n",file=outfile,append=T)
		cat("----------------------------------\n",file=outfile,append=T)
		cat("-tips ....",all(tipLabels(test1) == tipLabels(test1.remix)),"\n")
		cat("----------------------------------\n",file=outfile,append=T)		
		cat("weights: ",weight(test1),"\n",file=outfile,append=T)
		cat("weights: ",weight(test1.remix),"\n",file=outfile,append=T)	
		
	}
	
	# EXECUTE
	
}

cat("\n\nStarting execution tests....",file=outfile,append=T)

for(ff in nexfiles)
{
	cat(ff,"\n",file=outfile,append=T)
	read.brownie(paste(ff,".tmp",sep=""))
}

