# To do contour plot of brlen in R: #

  1. copy the table of brlen from besttrees.tre. The first col is the ∆lnL, the others are the brlen
  1. save this in a file (say, Untitled.txt)
  1. to do the contour plot on log of the brlen in R:
```
a<-read.table("Untitled.txt");
library(akima);
data.interp2 <- interp(log(a[,2]), log(a[,3]), a[,1]); #here, take out log if you don't want to plot log of brlen
filled.contour(data.interp2,levels=c(0,0.000000000000001,2,4,6,8,10,12),col=gray(c(0:8)/8)); #the odd 0.000..1 level is so you can see the area where the lnL of the brlen are all equal
```