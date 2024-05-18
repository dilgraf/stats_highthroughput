library(devtools)
install_github("genomicsclass/GSE5859Subset")
#High-throughput technologies measure thousands of _features_,
#_samples_ - experimental units
#So a high-throughput experiment is usually defined by three tables: 
#one with the high-throughput measurements and 
#two tables with information about the columns and rows of this first table.
# transposed n*m matrix, where n - the number of units and m - the number of features

library(GSE5859Subset)
data(GSE5859Subset)
head(geneExpression) #features(genes [ids provided by manufacturer])*units(samples)
dim(geneExpression)#8793*24

head(sampleInfo) #24*4
identical(colnames(geneExpression), sampleInfo$filename) #TRUE

#How many samples were processed on 2005-06-27?
sampleInfo[sampleInfo$date == "2005-06-27",]

#How many of the genes are on chromosome Y?
geneAnnotation[(geneAnnotation$CHR == "chrY")&(!is.na(geneAnnotation$CHR)),]

#log expression value for gene ARPC1A that was measured on 2005-06-10?
geneExpression[which(geneAnnotation$SYMBOL=="ARPC1A"),
               which(sampleInfo$date == "2005-06-10")]

#median value of each column
median(apply(geneExpression, 2, median))

#smallest p-value among t-tests for cases (1) and controls (0):
min(apply(geneExpression, 1, function(x) 
  t.test(x[which(sampleInfo$group==1)], x[which(sampleInfo$group==0)])$p.value))


# Looking for differentially expressed genes (cases' average is not the same as for controls)
g <- sampleInfo$group
#simplifying our experiment, we'll look only into 1 feature (gene):
e <- geneExpression[25,]
#before applying ttest, performing normality check:
par(mfrow=c(1,2))
qqnorm(e[g==1])
qqline(e[g==1])
qqnorm(e[g==0])
qqline(e[g==0])

#t-test:
t.test(e[g==1], e[g==0])
#p>0.05 - the difference not significant
mytest <- function(x){
  t.test(x[g==1], x[g==0], var.equal=TRUE)$p.value
}
pvals <- apply(geneExpression, 1, mytest)
sum(pvals <= 0.05)#1383

#creating a matrix using Monte Carlo simulation
m <- nrow(geneExpression)
n <- ncol(geneExpression)
randomData <- matrix(rnorm(n*m),m,n)
#applying t-test to each row of this new matrix for which we
#know the null hypothesis is true for every single feature.
nullpval <- apply(randomData, 1, mytest)
sum(nullpval<= 0.05)#452 significant genes, where we should've gotten 0

#Now we will run these 20 experiments 1,000 times and each time 
#save the number of p-values that are less than or equal to 0.05
set.seed(100)
simulate <- function() {
  pvals <- replicate(20,{
  control = rnorm(10,30,2)
  treatment = rnorm(10,30,2)
  t.test(treatment,control)$p.value
  })
  sum(pvals<=0.05)
}
p_vals_below <- replicate(1000, simulate())
mean(p_vals_below)

#for what proportion of the 1,000 replicates do we reject the null hypothesis 
#at least once (more than 0 false positives)
sum(p_vals_below>0)/length(p_vals_below)