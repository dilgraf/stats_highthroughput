filename <- 'femaleControlsPopulation.csv'
#controls:
population <- unlist(read.csv(filename))

#Goal: to simulate an example where we try 10,000 different diets to see
#if they have an effect on weight.
set.seed(1)
#We're just going to generate samples (n=12) as if H0 is true for all diets
N <- 12
m <- 10000
pvals <- replicate(m,{
  control = sample(population, N)
  treatment = sample(population, N)
  t.test(treatment,control)$p.value
})
#R:
sum(pvals < 0.05)
#in this case R=V (H0 is true)

###Alternative is true for 10% of diets
alpha <- 0.05
p0 <- 0.90 #10% of diets work, 90% - don't
m0 <- m*p0
m1 <- m-m0
nullHypothesis <- c(rep(TRUE, m0), rep(FALSE, m1))
delta <- 3 #difference between two diets by 3 grams

set.seed(1)
calls <- sapply(1:m, function(i){
  control = sample(population, N)
  treatment = sample(population, N)
  #when null hypothesis is false, we add delta:
  if(!nullHypothesis[i]) treatment <- treatment+delta
  ifelse(t.test(treatment,control)$p.value < alpha,
         "Called significant",
         "Not called significant")
})

##We usually don't know when the H0 is true, but here we do.
null_hypothesis <- factor(nullHypothesis, levels = c("TRUE", "FALSE"))
table(null_hypothesis, calls)

#we can run the experiment multiple times and the number will change
#number of simulations:
B <- 10
VandS <- replicate(B,{
  calls <- sapply(1:m, function(i){
    control = sample(population, N)
    treatment = sample(population, N)
    #when null hypothesis is false, we add delta:
    if(!nullHypothesis[i]) treatment <- treatment+delta
    t.test(treatment,control)$p.value < alpha
    })
  cat("V =", sum(nullHypothesis & calls), 
      "S =", sum(!nullHypothesis & calls), "\n")
  c(sum(nullHypothesis & calls), sum(!nullHypothesis & calls))
})
#This procedure has a very high familywise error rate.

#vectorizing code to make it faster 
#(In R, operations based on matrices are typically much faster 
#than operations performed within loops or sapply())
library(genefilter)
set.seed(1)
g <- factor(c(rep(0,N), rep(1,N)))
B <- 10
system.time(
  VandS <- replicate(B, {
    #matrix with control data (rows - tests, columns - mice)
    controls <- matrix(sample(population, N*m, replace = TRUE), nrow=m)
    treatments <- matrix(sample(population, N*m, replace = TRUE), nrow=m)
    #add effect (+3g) to 10% of them
    treatments[which(!nullHypothesis),] <- treatments[which(!nullHypothesis),]+delta
    #combine to form one matrix
    dat <- cbind(controls, treatments)
    calls <- rowttests(dat, g)$p.value < alpha
    c(sum(nullHypothesis&calls), sum(!nullHypothesis&calls))
  })
)

##Bonferroni
controls <- matrix(sample(population, N*m, replace = TRUE), nrow=m)
treatments <- matrix(sample(population, N*m, replace = TRUE), nrow=m)
treatments[which(!nullHypothesis),] <- treatments[which(!nullHypothesis),]+delta
#combine to form one matrix
dat <- cbind(controls, treatments)
g <- factor(c(rep(0,N), rep(1,N)))
pvals <- rowttests(dat, g)$p.value
sum(pvals<= alpha/m) #only 2 cases are significant
#very high false negative rate

##FDR (false discovery rate)
calls <- pvals<=alpha
R <- sum(calls)
V <- sum(nullHypothesis&calls)
Q <- ifelse(R>0, V/R, 0) 
B <- 1000
Qs <- replicate(B,{
  controls <- matrix(sample(population, N*m, replace = TRUE), nrow=m)
  treatments <- matrix(sample(population, N*m, replace = TRUE), nrow=m)
  treatments[which(!nullHypothesis),] <- treatments[which(!nullHypothesis),]+delta
  #combine to form one matrix
  dat <- cbind(controls, treatments)
  calls <- rowttests(dat, g)$p.value < alpha
  R=sum(calls)
  Q=ifelse(R>0, sum(nullHypothesis&calls)/R, 0)
  return(Q)
})
hist(Qs)
FDR <- mean(Qs)
print(FDR)
#the FDR is high

fdr <- p.adjust(pvals, method = "fdr")
par(mfrow=c(1,1))
plot(pvals, fdr, log="xy")
#Benjamini-Hochberg procedure
alpha <- 0.05
B <- 100
res <- replicate(B,{
  controls <- matrix(sample(population, N*m, replace = TRUE), nrow=m)
  treatments <- matrix(sample(population, N*m, replace = TRUE), nrow=m)
  treatments[which(!nullHypothesis),] <- treatments[which(!nullHypothesis),]+delta
  #combine to form one matrix
  dat <- cbind(controls, treatments)
  pvals <- rowttests(dat, g)$p.value
  calls <- p.adjust(pvals, method = "fdr") < alpha
  R <- sum(calls)
  Q <- ifelse(R>0, sum(nullHypothesis&calls)/R, 0)
  return(c(R,Q))
})
Qs <- res[2,]
hist(Qs)
FDR = mean(Qs)
print(FDR)

###FDR exercises
library(devtools)
library(GSE5859Subset)
data(GSE5859Subset)
library(genefilter)
##1
pvals <- rowttests(geneExpression, factor(sampleInfo$group))$p.value
sum(pvals< 0.05)
##2 - Bonferroni
sum(pvals<= 0.05/length(pvals))
##3 - FDR
sum(p.adjust(pvals, method = "fdr")<0.05)
#Note that we are controlling two very different error rates. 
#Here we are saying that we think this list of 13 genes has about 5% FP. 
#The Bonferroni procedure gave us a list ot 10 genes for which we were quite certain had no FP. 
#Note again that we have not discussed FN.

##4 - qvalue
library(qvalue)
sum(qvalue(pvals)$qvalues<0.05)
#qvalue, estimates FDR differently and is less conservative. 
#Remember that the theory provides bounds for FDR: it guarantees FDR will be less than 0.05. 
#If qvalue does in fact estimate pi0 well then it will provide a list with FDR closer to 0.05.

##5 - pi0
qvalue(pvals)$pi0

##6
plot(qvalue(pvals)$qvalue/p.adjust(pvals,method="fdr"))
abline(h=qvalue(pvals)$pi0,col=2)

hist(pvals,breaks=seq(0,1,len=21))
expectedfreq <- length(pvals)/20 #per bin
abline(h=expectedfreq*qvalue(pvals)$pi0,col=2,lty=2)

##7
n <- 24
m <- 8793
mat <- matrix(rnorm(n*m),m,n)
delta <- 2
#H1 is true for 500 genes
positives <- 500
mat[1:positives,1:(n/2)] <- mat[1:positives,1:(n/2)]+delta

set.seed(1)
#FP using Bonferroni
res <- replicate(1000,{
  mat <- matrix(rnorm(n*m),m,n)
  mat[1:positives,1:(n/2)] <- mat[1:positives,1:(n/2)]+delta
  pvals <- rowttests(mat, factor(sampleInfo$group))$p.value
  #only including p-values of genes for which Ho is true
  FP <- sum(pvals[-(1:positives)] <= 0.05/m)
})
mean(res/8293)
#Bonferroni controls FWER to be 0.05 not the FDR -> FDR is extremely low

##8
#FN using Bonferroni
res <- replicate(1000,{
  mat <- matrix(rnorm(n*m),m,n)
  mat[1:positives,1:(n/2)] <- mat[1:positives,1:(n/2)]+delta
  pvals <- rowttests(mat, factor(sampleInfo$group))$p.value
  FN <- sum(pvals[1:positives] > 0.05/m)
})
mean(res/500)
#low FDR -> FN is high

##9
#FP using qvalues from padjust
res <- replicate(1000,{
  mat <- matrix(rnorm(n*m),m,n)
  mat[1:positives,1:(n/2)] <- mat[1:positives,1:(n/2)]+delta
  pvals <- rowttests(mat, factor(sampleInfo$group))$p.value
  padj <- p.adjust(pvals, method = "fdr")
  FP <- sum(padj[-(1:positives)] <= 0.05)
})
mean(res/8293)
#the Benjamini–Hochberg procedure gives us a bound. 
#The larger m1, the more conservative this approximation will be.

##10
#FN using padjusted qvalues
res <- replicate(1000,{
  mat <- matrix(rnorm(n*m),m,n)
  mat[1:positives,1:(n/2)] <- mat[1:positives,1:(n/2)]+delta
  pvals <- rowttests(mat, factor(sampleInfo$group))$p.value
  padj <- p.adjust(pvals, method = "fdr")
  FP <- sum(padj[1:positives] > 0.05)
})
mean(res/500)
#the potential advantage of FDR over FWER, the FN is much reduced

##11
#FP using qvalues from qvalues
res <- replicate(1000,{
  mat <- matrix(rnorm(n*m),m,n)
  mat[1:positives,1:(n/2)] <- mat[1:positives,1:(n/2)]+delta
  pvals <- rowttests(mat, factor(sampleInfo$group))$p.value
  qvals <- qvalue(pvals)$qvalues
  FP <- sum(qvals[-(1:positives)] <= 0.05)
})
mean(res/8293)
#by estimating pi0 this approach gets closer to the targeted FDR of 0.05

##12
#FN using qvalues from qvalue
res <- replicate(1000,{
  mat <- matrix(rnorm(n*m),m,n)
  mat[1:positives,1:(n/2)] <- mat[1:positives,1:(n/2)]+delta
  pvals <- rowttests(mat, factor(sampleInfo$group))$p.value
  qvals <- qvalue(pvals)$qvalues
  FP <- sum(qvals[1:positives] > 0.05)
})
mean(res/500)
#by creating a list of an FDR closer to 0.05 we are less conservative and 
#thus decrease the false negative rate further