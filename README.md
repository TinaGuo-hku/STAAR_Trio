# STAAR_Trio

## Setup
Please use the following command in R to install the package:
```
install.packages("devtools") 
library(devtools)
install_github("TinaGuo-hku/STAAR_Trio")
```
We generated 10 functional annotations and calculated the Phred score of the 10 annotations

```
numVar <- dim(snploc)[1]
Z1 <- rnorm(numVar)
Z2 <- rnorm(numVar)
Z3 <- rnorm(numVar)
Z4 <- rnorm(numVar)
Z5 <- rnorm(numVar)
Z6 <- rnorm(numVar)
Z7 <- rnorm(numVar)
Z8 <- rnorm(numVar)
Z9 <- rnorm(numVar)
Z10 <- rnorm(numVar)
Z <- cbind(Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9,Z10)

rank <- colRanks(Z,preserveShape = TRUE)
PHRED <- -10*log10(1-rank/dim(rank)[1])

```
Then randomly selected a 5kb signal region from the 1-Mb region.

```
# Length of the sequence
maxLength <- 1000000
# Length of the signal region
sigLength <- 5000
# Location of the signal region
startloc <- sample(1:(maxLength-sigLength+1),1)
endloc <- startloc + sigLength - 1
# Extract the signal region from snploc
snplist <- which(snploc$CHROM_POS >= startloc & snploc$CHROM_POS <= endloc)
snpRegion <- snploc[snplist,]
# number of variants in the signal region
numSNPs <- dim(snpRegion)[1]
# genotype of the signal region
Geno <- genotype[,snplist]
```
We first randomly select 5 annotations to determine the probability of variants being causal in the region through a logistic model. We then randomly selected casual variants based on these probabilities. 

```
# Simulate causal variants from this region

b0 <- rje::logit(0.18)
b <- rep(log(5),10)

causalprob <- apply(Z[snplist,],1,function(z){
  ind <- sample(1:10,5)
  rje::expit(b0 + b[ind] %*% z[ind])
})
isCausal <- rbinom(numSNPs,1,causalprob)

```
Calculate the p values of the simulated data using STAAR_Trio


Here beta=-0.13*log10(MAF). The directions were all positive. The phenotype of samples was generated through a linear model.

```
# effects size of covariates 
alpha0 <- 0.5
# effect size of causal variants
c0 <- 0.13
beta <- -c0 * log10(maf[snplist][which(isCausal==1)])
# generate phenotype data
Y <- alpha0 + as.matrix(genotype)[,snplist][,which(isCausal==1)] %*% beta + eps
# offspring phenotype data
Y_off <- Y[seq(3,15000,3)]

```
We applied STAAR_Trio to analyze the simulated data set.

```
# calcualte the residuals 
null_model <- lm(Y~1)
staar_trio(trio[id_trio,snplist],maf[snplist],pos[snplist],PHRED[snplist,],adjust_for_covariates=TRUE,y_res=null_model$residuals[seq(3,15000,3)])
```




