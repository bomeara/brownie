rm(list=ls())
# Test setup:
require(RBrownie)
require(grid)


cmptrees <- function(t1,t2)
{
	retval = TRUE
	
	if(nNodes(t1) != nNodes(t2)) {
		warning("nNodes don't match")
		retval = FALSE
	}
	
	if(nTips(t1) != nTips(t2)) {
		warning("nTips don't match")
		retval = FALSE
	}
	
	if(nEdges(t1) != nEdges(t2)) {
		warning("nEdges don't match")
		retval = FALSE
	}
	
	if(xor(hasSingle(t2),hasSingle(t1))) {
		warning("hasSingle don't match")
		retval = FALSE
	}
	
	if(xor(isRooted(t2),isRooted(t1))) {
		warning("isRooted don't match")
		retval = FALSE
	}
			
	if(rootNode(t1) != rootNode(t2)) {
		warning("rootNode don't match")
		retval = FALSE
	}
	
	# edge lengths:
	if(any(is.na(edgeLength(t1))))
		t1@edge.length[is.na(edgeLength(t1))] <- 0.000
		
	if(any(is.na(edgeLength(t2))))
		t2@edge.length[is.na(edgeLength(t2))] <- 0.000

			
	if(any(sort(round(edgeLength(t1),2)) != sort(round(edgeLength(t2),2)))) {
		warning("rootNode don't match")
		retval = FALSE
	}	
		
	retval
}


cmpSubNodes<-function(t1,t2)
{
	retval=TRUE
	
	if(xor(hasSubNodes(t2),hasSubNodes(t1))) {
		warning("hasSubNodes don't match")
		retval = FALSE
	}

	for(ii in seq(ncol(sndata(t1))))
	{
		if(!all(sort(sndata(t1)[,ii]) == sort(sndata(t2)[,ii]))) {
			warning("sndata don't match",ii)
			retval = FALSE
		}
	}
	
	if(!all(round(sort(apply(snposition(t1),1,mean))) == round(sort(apply(snposition(t2),1,mean))))) {
		warning("snposition don't match")
		retval = FALSE
	}

	if(!all(sort(apply(snbranch(t1),1,sum)) == sort(apply(snbranch(t2),1,sum)))) {
		warning("snbranch don't match")
		retval = FALSE
	}

	retval
}

p1=proc.time()

## test block 1:
test1 = "(('Taxon1':{A,0.2; A,.5; C,0.1}, 'Taxon2':{T,0.2; C,0.15}):{C,0.5}, Taxon3:{C,0.4} ):{G,0};"

# is simmap tests
stopifnot(is.simmap(text=test1,vers=1.0))
stopifnot(!is.simmap(text=test1,vers=1.5))

stopifnot(is.simmap(text=test1,vers=1.0,quick=FALSE))
stopifnot(!is.simmap(text=test1,vers=1.5,quick=FALSE))

# (intitial) read and convert tests:
test1.tree = read.simmap(text=test1,vers=1.0)
test1.ext=as(test1.tree,"phylo4d_ext")
stopifnot(is(test1.ext, "phylo4d_ext"))
test1.ext=phyext(test1.tree)
stopifnot(is(test1.ext, "phylo4d_ext"))


# manipulate tests (add/remove/edit subnodes)

# read tests
test1.read1 = try(read.simmap(text=test1,vers=1.0),silent=T)
test1.read2 = try(read.simmap(text=test1,vers=1.1),silent=T)
test1.read3 = try(read.simmap(text=test1,vers=1.5),silent=T)
stopifnot(inherits(test1.read1,'phylo4d'))
stopifnot(is(test1.read2,'try-error'))
stopifnot(is(test1.read3,'try-error'))

# write tests:
test1.outted1 = capture.output(write.simmap(test1.ext,file=stdout(),vers=1.0))
test1.reread1 = read.simmap(text=test1.outted1,vers=1.0)
stopifnot(cmptrees(test1.reread1,test1.tree))
stopifnot(cmpSubNodes(phyext(test1.reread1),test1.ext))
test1.outted2 = capture.output(write.simmap(test1.ext,file=stdout(),vers=1.1))
test1.reread2 = read.simmap(text=test1.outted2,vers=1.1)
stopifnot(cmptrees(test1.reread2,test1.tree))
stopifnot(cmpSubNodes(phyext(test1.reread2),test1.ext))
test1.outted3 = capture.output(write.simmap(test1.ext,file=stdout(),vers=1.5))
test1.reread3 = read.simmap.new(text=test1.outted3)
stopifnot(cmptrees(test1.reread3,test1.ext))
stopifnot(cmpSubNodes(test1.reread3,test1.ext))



# plot tests
plot(test1.ext)
## end
p2 = proc.time() - p1
print(p2)


## test block 1.5:
data(geospiza)
test15.tree = geospiza
test15.ext=as(test15.tree,"phylo4d_ext")
stopifnot(is(test15.ext, "phylo4d_ext"))
test15.ext=phyext(test15.tree)
stopifnot(is(test15.ext, "phylo4d_ext"))

# read/write tests
test15.outted1 = capture.output(write.simmap(test15.ext,file=stdout(),vers=1.5))
test15.reread1 = read.simmap.new(text=test15.outted1)
stopifnot(cmptrees(test15.reread1,test15.tree))
stopifnot(cmpSubNodes(phyext(test15.reread1),test15.ext))

# manipulate subnodes
# (write 1.1 and 1.0 format, then remove data from test15.ext

## end


## test block 2:
test2 = "(((((32:{0,0.058907691963},(31:{0,0.022424212177:1,0.016767967599},(34:{1,0.036771632295},(30:{1,0.033928324112},33:{0,0.029725747655:1,0.004202576457}):{1,0.002843308183}):{1,0.002420547480}):{1,0.000275233450:0,0.019440278711}):{0,0.008552585573},29:{0,0.067460277504}):{0,0.040135207374},8:{0,0.091158418718:1,0.000880105051:0,0.015556961136}):{0,0.019704981181},(10:{0,0.112116752589},(2:{0,0.000772585273},3:{0,0.000772585273}):{0,0.014969340499:1,0.029030780722:0,0.067344046095}):{0,0.015183713480}):{0,0.023742583944},((1:{1,0.066050466671},4:{1,0.066050466671}):{1,0.039823924285},(9:{1,0.086482781129},(((5:{1,0.017325369522},6:{1,0.017325369522}):{1,0.006182126280},7:{1,0.023507495803}):{1,0.034790702136},(((26:{1,0.019800202825},14:{1,0.019800202825}):{1,0.016420404516},22:{1,0.036220607341}):{1,0.009769629593},(((23:{1,0.014983924906},(24:{1,0.004248281035},(18:{1,0.001668719645},21:{1,0.001668719645}):{1,0.002579561389}):{1,0.010735643875}):{1,0.009021350133},12:{1,0.024005275044}):{1,0.015472368742},(((((20:{1,0.019149815099},(15:{1,0.015578674920},11:{1,0.015578674920}):{1,0.003571140182}):{1,0.006407372276},25:{1,0.025557187375}):{1,0.001278231913},(13:{1,0.017974442632},16:{1,0.017974442632}):{1,0.008860976659}):{1,0.004074632240},(17:{1,0.027243985323},(28:{1,0.022983956116},19:{1,0.022983956116}):{1,0.004260029207}):{1,0.003666066205}):{1,0.005572724711},27:{1,0.036482776240}):{1,0.002994867546}):{1,0.006512593148}):{1,0.012307961005}):{1,0.028184583180}):{1,0.019391609822}):{1,0.019245403638:0,0.025923255434});"
test2.tree = read.simmap(text=test2)
## end


## test block 3:
test3 = "((38:[&map={0}]0.032983,(36:[&map={0}]0.028303,37:[&map={1,0.010533,0}]0.028303):[&map={0}]0.004680):[&map={0}]0.050180,(1:[&map={1,0.012333,0}]0.063472,(((7:[&map={1}]0.009876,8:[&map={1}]0.009876):[&map={1,0.022100,0}]0.029626,(6:[&map={1}]0.018863,((2:[&map={0}]0.004847,3:[&map={0}]0.004847):[&map={0,0.002049,1}]0.008509,(4:[&map={0,0.000652,1}]0.001499,5:[&map={1}]0.001499):[&map={1}]0.011857):[&map={1}]0.005507):[&map={1,0.005990,0}]0.020639):[&map={0}]0.000317,((34:[&map={0}]0.030583,35:[&map={1,0.005646,0}]0.030583):[&map={0}]0.005200,(((32:[&map={0}]0.012017,33:[&map={0}]0.012017):[&map={0}]0.023234,(31:[&map={0}]0.031836,((24:[&map={0}]0.028630,(22:[&map={0}]0.011680,23:[&map={1,0.004353,0}]0.011680):[&map={0,0.010175,1,0.005622,0}]0.016950):[&map={0}]0.002796,((25:[&map={0}]0.021272,26:[&map={0}]0.021272):[&map={0}]0.008918,(30:[&map={0}]0.026957,(27:[&map={1,0.013033,0}]0.026116,(28:[&map={0}]0.023870,29:[&map={0}]0.023870):[&map={0}]0.002246):[&map={0}]0.000841):[&map={0}]0.003233):[&map={0}]0.001236):[&map={0}]0.000410):[&map={0}]0.003416):[&map={0}]0.000360,((17:[&map={0}]0.019696,(18:[&map={0}]0.015763,(((19:[&map={1}]0.001483,39:[&map={1}]0.001483):[&map={1}]0.002346,40:[&map={1}]0.003829):[&map={1}]0.004309,(20:[&map={0,0.000209,1}]0.001346,21:[&map={1}]0.001346):[&map={1}]0.006792):[&map={1,0.005543,0}]0.007625):[&map={0}]0.003933):[&map={0}]0.015817,((9:[&map={1}]0.000782,10:[&map={1}]0.000782):[&map={1,0.024029,0}]0.032488,((15:[&map={0}]0.011159,16:[&map={0}]0.011159):[&map={0}]0.011243,(14:[&map={0}]0.021356,(13:[&map={0}]0.020872,(11:[&map={0}]0.004325,12:[&map={0}]0.004325):[&map={0}]0.016547):[&map={0}]0.000484):[&map={0}]0.001045):[&map={0}]0.010869):[&map={0}]0.002242):[&map={0}]0.000098):[&map={0}]0.000172):[&map={0}]0.004037):[&map={0}]0.023653):[&map={0}]0.019690);"
test3.tree = read.simmap.new(text=test3)


## end test block 3


## test block 4:
test4.fname = "simmap50trees.txt"
test4.trees = read.nexus.simmap(test4.fname)
stopifnot(length(test4.trees) == 50)
write.nexus.simmap(test4.trees,file="simmap50trees_OUT.txt",vers=1.5)
test4.reread = read.nexus.simmap("simmap50trees_OUT.txt")

for(ii in seq(length(test4.trees)))
	stopifnot(cmptrees(test4.reread[[ii]],test4.trees[[ii]]))

##


## test block 5:
test5.fname = "CCHF_S.trees"
test5.trees = read.nexus.simmap(test5.fname)
stopifnot(length(test5.trees) == 3001)
test5.trees = phyext(test5.trees)
write.nexus.simmap(test5.trees,file="CCHF_S_OUT.trees",vers=1.5)
test5.reread = read.nexus.simmap("CCHF_S_OUT.trees")

for(ii in seq(length(test5.trees)))
	stopifnot(cmptrees(test5.reread[[ii]],test5.trees[[ii]]))


##

