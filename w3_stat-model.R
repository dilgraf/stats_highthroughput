##Poisson distribution
N=10000 #number of genes
lambdas=2^seq(1,16,len=N) #these are the true abundances of genes
#generate two replicates
#H0 is true for all genes
y=rpois(N, lambdas)
x=rpois(N, lambdas)
#if it's appropriate to consider the log ratio as a summary statistic?
#Log ratios are commonly used in biology to summarize differences
ind = which(y>0 & x>0)    # make sure no 0s due to ratio and log

#plot the log ratio against the lambdas
plot(log2(lambdas), log2(y/x))
#for low values of lambda, you can get very large log ratios

##real gene expression data
library(parathyroidSE)
data(parathyroidGenesSE)
se = parathyroidGenesSE
x = assay(se)[,23]
y = assay(se)[,24]
ind = which(y>0 & x>0)
#since we don't have lambda, we'll use the average
plot((log2(x)+log2(x))/2, log2(y/x))

##Exercises
#1 - The probability of conceiving a girl is 0.49. 
#What is the probability that a family with 4 children has 2 girls and 2 boys
dbinom(2, 4, 0.49)

#2 - probability that a family with 10 children has 4 girls and 6 boys
dbinom(4, 10, 0.49)

#3 - what is the probability that the GC-content (proportion of Gs or Cs) is strictly above 0.5 in this interval
1-pbinom(10, 20, 0.4)

#4 - what is the probability that at least one winning tickets is sold?
1-dbinom(0, 189000000, 1/175223510)

#5 - probability that two or more winning tickets are sold
1-pbinom(1, 189000000, 1/175223510)

#6 - what is the exact probability that the GC-content (proportion of Gs of Cs) 
#is greater than 0.35 and smaller or equal to 0.45 in this interval
pbinom(20*0.45, 20, 0.4)-pbinom(20*0.35, 20, 0.4)

#7 - what is the normal approximation to the probability?
pnorm((20*0.45-20*0.4)/sqrt(20*0.4*0.6))-pnorm((20*0.35-20*0.4)/sqrt(20*0.4*0.6))

#8
exact = pbinom(1000*0.45, 1000, 0.4)-pbinom(1000*0.350, 1000, 0.4)
approxim = pnorm((1000*0.45-1000*0.4)/sqrt(1000*0.4*0.6))-pnorm((1000*0.35-1000*0.4)/sqrt(1000*0.4*0.6))
abs(exact-approxim)

#9
Ns <- c(5,10,30,100)
ps <- c(0.01,0.10,0.5,0.9,0.99)
par(mfrow=c(4,5))
for (N in Ns){
  k <- 1:(N-1)
  for (p in ps){
    exact = dbinom(k,N,p)
    a <- (k+0.5 - N*p)/sqrt(N*p*(1-p))
    b <- (k-0.5 - N*p)/sqrt(N*p*(1-p))
    approx = pnorm(a) - pnorm(b)
    LIM <- range(c(approx,exact))
    plot(exact,approx, xlim=LIM,ylim=LIM, col=1,pch=16)
    abline(0,1)
  }
}

#10 - Poisson
N <- 189000000
p <- 1/175223510
dbinom(2,N,p)
a <- (2+0.5 - N*p)/sqrt(N*p*(1-p))
b <- (2-0.5 - N*p)/sqrt(N*p*(1-p))
pnorm(a) - pnorm(b)
dpois(2,N*p) #practically the same because  is very very large and  is not 0
1-ppois(1,N,p)

##MLE - Maximum Likelihood Estimate
library(devtools)
#install_github("genomicsclass/dagdata")
library(dagdata)
data(hcmv)
breaks=seq(0,4000*round(max(locations)/4000),4000)
#cut - converts numeric to factor
tmp=cut(locations,breaks)
counts=as.numeric(table(tmp))
l<-function(lambda) {
  ls <- dpois(counts,lambda,log=TRUE)
  return(sum(ls))
}
lambdas <- seq(3,7,len=100)
ls <- exp(sapply(lambdas,l))
plot(lambdas,ls,type="l")
mle=optimize(l,c(0,10),maximum=TRUE)
abline(v=mle$maximum)
print(c(mle$maximum, mean(counts)))

##1
#locations of palindromes
plot(locations,rep(1,length(locations)))
#These palindromes are quite rare, p is very small
#If we break the genome into bins of 4000 basepairs, 
#then we have Np not so small and we might be able to use Poisson
breaks=seq(0,4000*round(max(locations)/4000),4000)
tmp=cut(locations,breaks)
#number of palindromes in each bin
counts=as.numeric(table(tmp))
hist(counts)

#for lambda=4:
probs <- dpois(counts,4)
likelihood <- prod(probs)
likelihood
logprobs <- dpois(counts,4,log=TRUE)
loglikelihood <- sum(logprobs)
loglikelihood

lambdas = seq(0,15,len=300)
l<-function(lambda, x) {
  sum(ls <- dpois(x,lambda,log=TRUE))
}
ls <- sapply(lambdas,function(lambda) l(lambda,counts))
plot(lambdas,ls)
mle=lambdas[which.max(ls)]
print(mle)

##2
binLocation=(breaks[-1]+breaks[-length(breaks)])/2
plot(binLocation,counts,type="l",xlab=)
binLocation[which.max(counts)]

##4 - probability of seeing a count of 14 or more
lambda = mean(counts[ - which.max(counts) ])
1-ppois(13, lambda)

##6 - bonferroni
0.05/57

##7
#x-axis values
ps <- (seq(along=counts) - 0.5)/length(counts)
poisq <- qpois(ps,lambda)
qqplot(poisq, counts)
abline(0,1)

##Models for variance
#install_github("genomicsclass/tissuesGeneExpression")
library(tissuesGeneExpression)
data("tissuesGeneExpression")
library(genefilter)
y = e[,which(tissue=="endometrium")]

##1
s2 <- rowVars(y)
qqnorm(s2)
qqline(s2)
qqnorm(sqrt(s2))
qqline(sqrt(s2))

##2
library("limma")
fitFDist(s2, 14)$scale

##3
ps <- (seq(along=s2) - 0.5)/length(s2)
#quantiles predicted by the F-distribution
qfs <- qf(ps,14, fitFDist(s2, 14)$df2)*fitFDist(s2, 14)$scale
lim <- sqrt(range(c(qfs, s2)))

qqplot(sqrt(qfs), sqrt(s2),ylim=lim, xlim=lim)
abline(0,1)

#excl upper 5%
k <- sqrt(quantile(s2,0.95))
qqplot(sqrt(qfs), sqrt(s2), ylim=c(0,k), xlim=c(0,k))
abline(0,1)
