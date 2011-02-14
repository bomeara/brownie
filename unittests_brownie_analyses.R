require(RBrownie)

# Discrete tests ("discrete")
junk = readBrownie("geospiza.nex")
junkrun=runDiscrete(junk,brfile="disc_junk.txt",
models=c("nonrev","nonrev","equal","rev","rev","nonrev"),
freqs=c("unif","equilib","empirical","unif","OPTIMIZE","empirical"),
reconstruct=T)
#junkrun
#plot(junkrun$trees[[1]])

# Noncensored ("cont")
junk = readBrownie("parrot.nex")
df1=runNonCensored(junk,brfile="cont_junk_all.txt",models=brownie.models.continuous()[1:5],treeloop=F,charloop=F)

#df2=runNonCensored(junk,brfile="cont_junk_all_mixed.txt",models=brownie.models.continuous()[sample(1:5)],treeloop=F,charloop=F)
writeBrownie(junk,file="parrot.nex.tmp")
junk2 = readBrownie("parrot.nex.tmp")
df2=runNonCensored(junk2,brfile="cont_junk_all.txt",models=brownie.models.continuous()[1:5],treeloop=F,charloop=F)


#summaryCont(junkrun)
#junkrun=runNonCensored(junk,brfile="cont_junk_noloop.txt",models=brownie.models.continuous()[1:5],treeloop=F,charloop=F)
#summaryCont(junkrun)
junkrun=runNonCensored(junk,brfile="cont_junk_noloop1.txt",models=brownie.models.continuous()[1],treeloop=F,charloop=F)
stopifnot(nrow(junkrun)!=0)
junkrun=runNonCensored(junk,brfile="cont_junk_noloop2.txt",models=brownie.models.continuous()[2],treeloop=F,charloop=F)
stopifnot(nrow(junkrun)!=0)
junkrun=runNonCensored(junk,brfile="cont_junk_noloop3.txt",models=brownie.models.continuous()[3],treeloop=F,charloop=F)
stopifnot(nrow(junkrun)!=0)
junkrun=runNonCensored(junk,brfile="cont_junk_noloop4.txt",models=brownie.models.continuous()[4],treeloop=F,charloop=F)
stopifnot(nrow(junkrun)!=0)
junkrun=runNonCensored(junk,brfile="cont_junk_noloop5.txt",models=brownie.models.continuous()[5],treeloop=F,charloop=F)
#stopifnot(nrow(junkrun)!=0) # this will be 0
junkrun=runNonCensored(junk,brfile="cont_junk_noloop6.txt",models=brownie.models.continuous()[1],treeloop=T,charloop=T)
stopifnot(nrow(junkrun)!=0)
junkrun=runNonCensored(junk,brfile="cont_junk_noloop7.txt",models=brownie.models.continuous()[2],treeloop=T,charloop=T)
stopifnot(nrow(junkrun)!=0)
junkrun=runNonCensored(junk,brfile="cont_junk_noloop8.txt",models=brownie.models.continuous()[3],treeloop=T,charloop=T)
stopifnot(nrow(junkrun)!=0)
junkrun=runNonCensored(junk,brfile="cont_junk_noloop9.txt",models=brownie.models.continuous()[4],treeloop=T,charloop=T)
stopifnot(nrow(junkrun)!=0)
junkrun=runNonCensored(junk,brfile="cont_junk_noloop10.txt",models=brownie.models.continuous()[5],treeloop=T,charloop=T)

junkrun=runNonCensored(junk,brfile="cont_junk_noloop_mixed.txt",models=brownie.models.continuous()[c(5,2,4,3,1)],treeloop=F,charloop=F)
summaryCont(junkrun) # TODO: fix!

# Censored ("ratetest")
junk = readBrownie("parrot.nex")
junkrun = runCensored(junk,brfile="ratetest_junk.txt",taxset="intrajoint",reps=1000,charloop=T)
summaryRatetest(junkrun)

