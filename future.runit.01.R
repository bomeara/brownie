require(RBrownie)
setwd("C:/Users/conrad/brownie")

# TEST 1
junk <- readBrownie("geospiza.nex")
datatypes(junk) == discData()

# TEST 2
junk2 <- readBrownie("parrot.nex")
datatypes(junk2[[1]]) == c(rep(contData(),2),rep(taxaData(),2))
taxasets(junk2[[1]]) <- sample(tipLabels(junk2[[1]]),10,replace=F)
datatypes(junk2[[1]]) == c(rep(contData(),2),rep(taxaData(),3))

# TEST 3
junk3 <- as(rtree(15),'brownie')
length(datatypes(junk3)) == 0
taxasets(junk3) <- sample(tipLabels(junk3),10,replace=F)
datatypes(junk3) == taxaData()
junkdata = data.frame(sample(c("A","b","D"),nTips(junk3),replace=T))
rownames(junkdata) <- tipLabels(junk3)
junk3 = addData(junk3,tip.data=junkdata)
datatypes(junk3) == c(taxaData(),discData())

# TEST 4/5
junk4 = as(rtree(15),'phylo4d')
!hasTipData(junk4)
taxasets(junk4) <- sample(tipLabels(junk4),10,replace=F)
taxasets(junk4) <- sample(tipLabels(junk4),10,replace=F)
hasTipData(junk4)
brownie(junk4)->junk5
datatypes(junk5) == rep(taxaData(),2)





