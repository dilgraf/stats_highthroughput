library(dplyr)
players <- read.csv("Batting.csv")
##2 - What is the average of these batting averages?
mean((filter(players,(yearID==2010) | (yearID==2011) | (yearID==2012)) %>% mutate(AVG=H/AB) %>% filter(AB>=500) %>% select(AVG))$AVG)

##3 - standard deviation
sd((filter(players,(yearID==2010) | (yearID==2011) | (yearID==2012)) %>% mutate(AVG=H/AB) %>% filter(AB>=500) %>% select(AVG))$AVG)

##5
sqrt((0.45*(1-0.45))/20)
##6
B<-0.11**2/(0.11**2+0.027**2)
0.275+(1-B)*(0.45-0.275)

library(Biobase)
library(SpikeInSubset)
data(rma95)
#extract the data matrix
y <- exprs(rma95)
#an experiment in which RNA was obtained from the same background pool 
#to create six replicate samples. Then RNA from 16 genes were artificially added 
#in different quantities to each sample
#These quantities (in picoMolars) and gene IDs are stored here
pData(rma95)

g <- factor(rep(0:1,each=3))
#create an index of which rows are associated with the artificially added genes
spike <- rownames(y) %in% colnames(pData(rma95))

##1 - What proportion of genes with a p-value < 0.01 are FP
library(genefilter)
tt <- rowttests(y, g)
mean(!spike[tt$p.value<0.01])

#volcano plot
smallp <- with(tt, p.value < .01)
cols <- ifelse(spike,"dodgerblue",ifelse(smallp,"red","black"))
with(tt, plot(-dm, -log10(p.value), cex=.8, pch=16, 
              xlim=c(-1,1), ylim=c(0,4.5),
              xlab="difference in means",
              col=cols))
abline(h=2,v=c(-.2,.2), lty=2)


##3
library(limma)
fit <- lmFit(y, design=model.matrix(~ g))
colnames(coef(fit))
fit <- eBayes(fit)
sampleSD = fit$sigma
posteriorSD = sqrt(fit$s2.post)
LIM = range( c(posteriorSD,sampleSD))
plot(sampleSD, posteriorSD,ylim=LIM,xlim=LIM)
abline(0,1)
abline(v=sqrt(fit$s2.prior))

##4
pvals = fit$p.value[,2]
mean(!spike[pvals<0.01])

with(tt, plot(fit$coef[,2], -log10(pvals), cex=.8, pch=16, 
                             xlim=c(-1,1), ylim=c(0,4.5),
                             xlab="difference in means",
                             col=cols))
abline(h=2,v=c(-.2,.2), lty=2)

#EDA for high-throughput experiments
library(GSE5859Subset)
data(GSE5859Subset)
g <- factor(sampleInfo$group)
results <- rowttests(geneExpression, g)

#simulation of p-values for which we know that H0 is true
m <- nrow(geneExpression)
n <- ncol(geneExpression)
randomData <- matrix(rnorm(n*m),m,n)
nullpvals <- rowttests(randomData, g)$p.value

## p-value histograms
par(mfrow=c(1,2))
hist(nullpvals) #uniform distribution
hist(results$p.value)

##volcano plot - shows us the effect size
#that we use in the t-test along with the p-values
#genes that have high significance are high on the plot
par(mfrow=c(1,1))
plot(results$dm,-log10(results$p.value),
     xlab="Effect size",ylab="- log (base 10) p-values")

##boxplots
library(Biobase)
#library(devtools)
#install_github("genomicsclass/GSE5859")
library(GSE5859)
data(GSE5859)
ge <- exprs(e)
ge[,49] <- ge[,49]/log2(exp(1)) ##immitate error
boxplot(ge,range=0,names=1:ncol(e),col=ifelse(1:ncol(ge)==49,1,2))
#we can simply show the boxplot summaries without the boxes:
qs <- t(apply(ge, 2, quantile, prob=c(0.05,0.25,0.5,0.75,0.95)))
matplot(qs,type="l", lty=1)

##MA plot
#to explore the relationship between two samples
x <- ge[,1]
y <- ge[,2]
par(mfrow=c(1,2))
plot(x,y)
#x-axis - average
plot((x+y)/2,x-y)
sd(y-x)
#standard deviation between these two measurements is 0.20 - not that small

##1
data(mas133)
e <- exprs(mas133)
par(mfrow=c(1,1))
plot(e[,1],e[,2],main=paste0("corr=",signif(cor(e[,1],e[,2]),3)),cex=0.5)
k <- 3000
b <- 1000 #a buffer
polygon(c(-b,k,k,-b),c(-b,-b,k,k),col="red",density=0,border="red")
#proportion of points inside
mean(e[,1]<k & e[,2]<k)

##2
plot(log2(e[,1]),log2(e[,2]),main=paste0("corr=",signif(cor(log2(e[,1]),log2(e[,2])),2)),cex=0.5)
k <- log2(3000)
b <- log2(0.5)
polygon(c(b,k,k,b),c(b,b,k,k),col="red",density=0,border="red")

##3
e <- log2(exprs(mas133))
plot((e[,1]+e[,2])/2,e[,2]-e[,1],cex=0.5)
sd(e[,2]-e[,1])
#if there is a mean shift the standard deviation will not summarize this. 
#We can instead consider the average distance:
sqrt(mean( (e[,2]-e[,1])^2))

##4
fold = abs(e[,2]-e[,1])
length(fold[fold>1])