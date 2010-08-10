require(RBrownie)

# Discrete tests ("discrete")
junk = readBrownie("geospiza.nex")
junkrun=runDiscrete(junk,brfile="disc_junk.txt",
models=c("nonrev","nonrev"),freqs=c("unif","equilib"),reconstruct=T)
junkrun
plot(junkrun$trees[[1]])

# Noncensored ("cont")
junk = readBrownie("parrot.nex")
junkrun=runNonCensored(junk,brfile="cont_junk.txt",models=brownie.models.continuous()[1:3],treeloop=T,charloop=T)
junkrun
sum.cont(junkrun)

# Censored ("ratetest")
junk = readBrownie("parrot.nex")
junkrun = runCensored(junk,taxset="intrajoint",reps=1000,charloop=T)
sum.ratetest(junkrun)

