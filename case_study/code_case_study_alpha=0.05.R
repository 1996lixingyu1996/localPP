# Data Analysis for BRAF-V600 study

# Load required packages
library(parallel)
library(doParallel)
library(R2jags)
library(foreach)
library(bhmbasket) 
library(basket) 
library(partitions)
library(plyr)
library(cluster)
library(coda)
library(rjags)
source("../functions/functions.R")
source("../functions/functions_sim.R")
source("../functions/functions_parallel.R")

nperclust <- 20000 # number of simulated trials per cluster 
detectCores()
nclust <- 5 # Total 100000 MC replicates
N <- rbind(19, 10, 26, 8, 14, 7
) #  total sample size for each indication without interim
B <- nrow(N) # total number of baskets
pnull <- rep(0.15, B) # null response rate for each indication
ptarget <- rep(0.35, B) # target response rate for each indication

## BOP2 for each indication with error rate sig.level=0.05
sig.level <- 0.05 # type I error
stopbounds <- NULL # since no interim
beta.a0 <- pnull # default beta prior 
beta.b0 <- 1-pnull # default beta prior 
ndigits = 3 ## number of digits for Q

## scenarios 
scenarios <- rbind( pnull ) # global null 

## generate data under global null to calibrate Q
seed <- 23235325
data.object<-generate.data(N, scenarios, ntrial = nperclust*nclust, seed = seed)

####################### Independent ########################
start_time <- Sys.time()
res.post <- post.infer.parallel(nclust = nclust, nperclust = nperclust, data.object = data.object, 
                                pnull = pnull, stopbounds =  stopbounds, beta.a0 = beta.a0, 
                                beta.b0 = beta.b0, seed = seed, ModelFit = "Independent")
Q <- get.Q.bwer(res.post, alpha = sig.level, digits = ndigits, Qclust = 1:B)
res <- get.weighted.power(res.post, Q = Q, alpha = sig.level)
end_time <- Sys.time()
(time.Independent <- end_time - start_time)
## Q
Q.Independent <- Q
## type I error
results.Independent <- list(res.post, res)
Qmat <- array(NA, dim = dim(res.post$postprob))
for(i in 1:B) Qmat[,,i] <- res$Q[i]
oc.Independent <- apply(res.post$postprob>Qmat, c(1,3), mean)
## posterior probabilities
fit <- Independent(nDat = c(19, 10, 26, 8, 14, 7), 
                   yDat = c(8, 0, 1, 1, 6, 2), 
                   be.a0 = beta.a0, be.b0 = beta.b0)
pp.Independent <- pbeta(0.15, shape1 = fit$a.post, 
                        shape2 = fit$b.post, lower.tail = FALSE)

####################### local PP ########################
a = 0.2
delta = 0.15
start_time <- Sys.time()
res.post <- post.infer.parallel(nclust = nclust, nperclust = nperclust, data.object = data.object, 
                                pnull = pnull, stopbounds =  stopbounds, beta.a0 = beta.a0, 
                                beta.b0 = beta.b0, seed = seed, ModelFit = "localPP", 
                                a = a, delta=delta)
Q <- get.Q.bwer(res.post, alpha = sig.level, digits = ndigits, Qclust = 1:B)
res <- get.weighted.power(res.post, Q = Q, alpha = sig.level)
end_time <- Sys.time()
(time.localPP <- end_time - start_time)
## Q
Q.localPP <- Q
## type I error
Qmat <- array(NA, dim = dim(res.post$postprob))
for(i in 1:B) Qmat[,,i] <- res$Q[i]
oc.localPP <- apply(res.post$postprob>Qmat, c(1,3), mean)
## posterior probabilities
fit <- localPP(nDat = c(19, 10, 26, 8, 14, 7), yDat = c(8, 0, 1, 1, 6, 2), 
               be.a0 = beta.a0, be.b0 = beta.b0, a = a, delta = delta)
pp.localPP <- pbeta(0.15, shape1 = fit$a.post, 
                        shape2 = fit$b.post, lower.tail = FALSE)
## similarity matrix
fit$sm # Table 7.3

#################################
### Table 4.2 Results on BRAF V600 trial data
################################
results <- matrix(NA, nrow = B, ncol = 6)
rownames(results) <- c("NSCLC", "CRC vemu", "CRC vemu+cetu", "Bile duct", "ECD or LCH", "ATC")
colnames(results) <- c("Q_ind", "Q_localPP", "TypeI_ind", "TypeI_localPP", "pp_ind", "pp_localPP")
results[,1] <- sprintf(fmt = '%#.3f', Q.Independent)
results[,2] <- sprintf(fmt = '%#.3f', Q.localPP)
results[,3] <- sprintf(fmt = '%#.3f', oc.Independent)
results[,4] <- sprintf(fmt = '%#.3f', oc.localPP)
results[,5] <- sprintf(fmt = '%#.3f', pp.Independent)
results[,6] <- sprintf(fmt = '%#.3f', pp.localPP)
results

save.image("case_study_alpha=0.05.RData")
