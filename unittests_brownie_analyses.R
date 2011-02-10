require(RBrownie)

# Discrete tests ("discrete")
junk = readBrownie("geospiza.nex")
junkrun=runDiscrete(junk,brfile="disc_junk.txt",
models=c("nonrev","nonrev","equal","rev","rev","nonrev"),
freqs=c("unif","equilib","empirical","unif","OPTIMIZE","empirical"),
reconstruct=T)
#junkrun
plot(junkrun$trees[[1]])

# Noncensored ("cont")
junk = readBrownie("parrot.nex")
junkrun=runNonCensored(junk,brfile="cont_junk_all.txt",models=brownie.models.continuous()[1:5],treeloop=T,charloop=T)
summaryCont(junkrun)
junkrun=runNonCensored(junk,brfile="cont_junk_all_mixed.txt",models=brownie.models.continuous()[sample(1:5)],treeloop=T,charloop=T)
summaryCont(junkrun)
junkrun=runNonCensored(junk,brfile="cont_junk_noloop.txt",models=brownie.models.continuous()[1:5],treeloop=F,charloop=F)
summaryCont(junkrun)
junkrun=runNonCensored(junk,brfile="cont_junk_noloop_mixed.txt",models=brownie.models.continuous()[c(5,2,4,3,1)],treeloop=F,charloop=F)
summaryCont(junkrun) # TODO: fix!

# Censored ("ratetest")
junk = readBrownie("parrot.nex")
junkrun = runCensored(junk,brfile="ratetest_junk.txt",taxset="intrajoint",reps=1000,charloop=T)
summaryRatetest(junkrun)

