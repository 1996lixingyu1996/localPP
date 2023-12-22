# Replicate the simulation results

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

nperclust <- 1000 # number of simulated trials per cluster 
detectCores()
nclust <- 5 # Total 5000 MC replicates
N <- rbind(c(20, 40),
           c(20, 40),
           c(20, 40),
           c(20, 40),
           c(20, 40)
) # interim sample size and total sample size for each indication
B <- nrow(N) # total number of baskets
pnull <- c(0.15, 0.15, 0.15, 0.15, 0.15) # null response rate for each indication
ptarget <- c(0.30, 0.30, 0.30, 0.30, 0.30) # target response rate for each indication

## BOP2 for each indication with error rate sig.level=0.1
sig.level <- 0.1 # type I error
stopbounds <- rbind(2, 
                    2, 
                    2, 
                    2, 
                    2) # obtained from BOP2 app
beta.a0 <- pnull # default beta prior 
beta.b0 <- 1-pnull # default beta prior 
ndigits = 3 ## number of digits for Q

## scenarios 
scenarios <- rbind( c(0.15, 0.15, 0.15, 0.15, 0.15),
                    c(0.15, 0.15, 0.15, 0.30, 0.30),
                    c(0.15, 0.30, 0.30, 0.30, 0.30),
                    c(0.15, 0.30, 0.30, 0.45, 0.45),
                    c(0.15, 0.45, 0.45, 0.45, 0.45),
                    c(0.30, 0.30, 0.30, 0.30, 0.30)
)

## generate data which will be analyzed by all methods
seed <- 887987
data.object<-generate.data(N, scenarios, ntrial = nperclust*nclust, seed = seed)

####################################################################################
## Selecting hyperparameter a for local PP
####################################################################################
aseq = seq(0, 1, 0.1)
delta <- 0.1
list.metrics0 <- list()
list.post <- list()
for(i in 1:length(aseq)){
  print(i)
  res.post <- post.infer.parallel(nclust = nclust, nperclust = nperclust, data.object = data.object, 
                                  pnull = pnull, stopbounds =  stopbounds, beta.a0 = pnull, 
                                  beta.b0 = 1-pnull, seed = seed, ModelFit = "localPP", 
                                  a = aseq[i], delta = delta)
  (Q0 <- get.Q.bwer(res.post, alpha = sig.level, digits = ndigits, Qclust = rep(1, B)))
  res0 <- get.weighted.power(res.post, Q = Q0, alpha = sig.level)
  list.metrics0[[i]] <- res0
  list.post[[i]] <- res.post
}
save.image("intermediate_results/localPP_select_a.RData")

####################################################################################
## Selecting hyperparameter tau for JSD
####################################################################################
aseq = seq(0, 1, 0.1) ## values for tau
list.metrics0 <- list()
list.post <- list()
for(i in 1:length(aseq)){
  print(i)
  res.post <- post.infer.parallel(nclust = nclust, nperclust = nperclust, data.object = data.object, 
                                  pnull = pnull, stopbounds =  stopbounds, beta.a0 = pnull, 
                                  beta.b0 = 1-pnull, seed = seed, ModelFit = "JSD", 
                                  tau = aseq[i])
  (Q0 <- get.Q.bwer(res.post, alpha = sig.level, digits = ndigits, Qclust = rep(1, B)))
  res0 <- get.weighted.power(res.post, Q = Q0, alpha = sig.level)
  list.metrics0[[i]] <- res0
  list.post[[i]] <- res.post
}
save.image("intermediate_results/JSD_select_tau.RData")


####################################################################################
## Simulation results for comparing the  following eight methods:
## Independent, local PP, JSD, EXNEX, BHM, BCHM, local MEM and MEM
####################################################################################

####################### Independent ########################
## start simulation
seed <- 887987
start_time <- Sys.time()
res.post <- post.infer.parallel(nclust = nclust, nperclust = nperclust, data.object = data.object, 
                                pnull = pnull, stopbounds =  stopbounds, beta.a0 = pnull, 
                                beta.b0 = 1-pnull, seed = seed, ModelFit = "Independent")
(Q <- get.Q.bwer(res.post, alpha = sig.level, digits = ndigits, Qclust = rep(1, B)))
res <- get.weighted.power(res.post, Q = Q, alpha = sig.level)
end_time <- Sys.time()
(time.Independent <- end_time - start_time)
## save results
results.Independent <- list(res.post, res)
Qmat <- array(NA, dim = dim(res.post$postprob))
for(i in 1:B) Qmat[,,i] <- res$Q[i]
(oc.Independent <- apply(res.post$postprob>Qmat, c(1,3), mean))
save.image("intermediate_results/AllMethods.RData")

####################### local PP ########################
## prior values
a = 0.2
delta = 0.1
## start simulation
seed <- 887987
start_time <- Sys.time()
res.post <- post.infer.parallel(nclust = nclust, nperclust = nperclust, data.object = data.object, 
                                pnull = pnull, stopbounds =  stopbounds, beta.a0 = pnull, 
                                beta.b0 = 1-pnull, seed = seed, ModelFit = "localPP", 
                                a = a, delta=delta)
(Q <- get.Q.bwer(res.post, alpha = sig.level, digits = ndigits, Qclust = rep(1, B)))
res <- get.weighted.power(res.post, Q = Q, alpha = sig.level)
end_time <- Sys.time()
(time.localPP <- end_time - start_time)
## save results
results.localPP <- list(res.post, res)
Qmat <- array(NA, dim = dim(res.post$postprob))
for(i in 1:B) Qmat[,,i] <- res$Q[i]
(oc.localPP <- apply(res.post$postprob>Qmat, c(1,3), mean))
save.image("intermediate_results/AllMethods.RData")

####################### JSD ########################
## prior values
tau = 0.5
## start simulation
seed <- 887987
start_time <- Sys.time()
res.post <- post.infer.parallel(nclust = nclust, nperclust = nperclust, data.object = data.object, 
                                pnull = pnull, stopbounds =  stopbounds, beta.a0 = pnull, 
                                beta.b0 = 1-pnull, seed = seed, ModelFit = "JSD", tau = tau)
(Q <- get.Q.bwer(res.post, alpha = sig.level, digits = ndigits, Qclust = rep(1, B)))
res <- get.weighted.power(res.post, Q = Q, alpha = sig.level)
end_time <- Sys.time()
(time.JSD <- end_time - start_time)
## save results
results.JSD <- list(res.post, res)
Qmat <- array(NA, dim = dim(res.post$postprob))
for(i in 1:B) Qmat[,,i] <- res$Q[i]
(oc.JSD <- apply(res.post$postprob>Qmat, c(1,3), mean))
save.image("intermediate_results/AllMethods.RData")

####################### lcoal MEM ########################
## prior values
delta = 2 
## start simulation
seed <- 887987
start_time <- Sys.time()
res.post <- post.infer.parallel(nclust = nclust, nperclust = nperclust, data.object = data.object, 
                                pnull = pnull, stopbounds =  stopbounds, beta.a0 = pnull, 
                                beta.b0 = 1-pnull, seed = seed, ModelFit = "localMEM", delta = delta)
(Q <- get.Q.bwer(res.post, alpha = sig.level, digits = ndigits, Qclust = rep(1, B)))
res <- get.weighted.power(res.post, Q = Q, alpha = sig.level)
end_time <- Sys.time()
(time.localMEM <- end_time - start_time)
## save results
results.localMEM <- list(res.post, res)
Qmat <- array(NA, dim = dim(res.post$postprob))
for(i in 1:B) Qmat[,,i] <- res$Q[i]
(oc.localMEM <- apply(res.post$postprob>Qmat, c(1,3), mean))
save.image("intermediate_results/AllMethods.RData")

####################### BHM uniform ########################
## prior values
u0 = 100
## start simulation
seed <- 887987
start_time <- Sys.time()
res.post <- post.infer.parallel(nclust = nclust, nperclust = nperclust, data.object = data.object, 
                                pnull = pnull, stopbounds =  stopbounds, beta.a0 = pnull, 
                                beta.b0 = 1-pnull, seed = seed, ModelFit = "BHMunif", u0 = u0)
(Q <- get.Q.bwer(res.post, alpha = sig.level, digits = ndigits, Qclust = rep(1, B)))
res <- get.weighted.power(res.post, Q = Q, alpha = sig.level)
end_time <- Sys.time()
(time.BHMunif <- end_time - start_time)
## save results
results.BHMunif <- list(res.post, res)
Qmat <- array(NA, dim = dim(res.post$postprob))
for(i in 1:B) Qmat[,,i] <- res$Q[i]
(oc.BHMunif <- apply(res.post$postprob>Qmat, c(1,3), mean))
save.image("intermediate_results/AllMethods.RData")

####################### BCHM ########################
## prior values
alpha = 1e-40 # precision for DP clustering
alpha1 = 50 # gamma(alpha1, beta1) for 1/variance
beta1 = 10 # gamma(alpha1, beta1) for 1/variance
mu0 = logit(pnull[1]) # hyperprior mean
tau2 = 0.01 # hyperprior Precision 
## start simulation
seed <- 887987
start_time <- Sys.time()
res.post <- post.infer.parallel(nclust = nclust, nperclust = nperclust, data.object = data.object, 
                                pnull = pnull, stopbounds =  stopbounds, beta.a0 = pnull, 
                                beta.b0 = 1-pnull, seed = seed, ModelFit = "BCHM", 
                                alpha = alpha, alpha1 = alpha1, beta1 = beta1, mu0 = mu0, tau2 = tau2)
(Q <- get.Q.bwer(res.post, alpha = sig.level, digits = ndigits, Qclust = rep(1, B)))
res <- get.weighted.power(res.post, Q = Q, alpha = sig.level)
end_time <- Sys.time()
(time.BCHM <- end_time - start_time)
## save results
results.BCHM <- list(res.post, res)
Qmat <- array(NA, dim = dim(res.post$postprob))
for(i in 1:B) Qmat[,,i] <- res$Q[i]
(oc.BCHM <- apply(res.post$postprob>Qmat, c(1,3), mean))
save.image("intermediate_results/AllMethods.RData")

####################### EXNEX ########################
## prior values
nex.w = 0.5
weight = c((1-nex.w)/2, (1-nex.w)/2, nex.w)
pguess = c(pnull[1], ptarget[1], (pnull[1] + ptarget[1])/2)
tau.HN.scale = rep(1, length(weight)-1)
Nmix <- length(weight)
Nexch <- Nmix - 1
mu.mean <- rep(NA, Nexch)
mu.prec <- rep(NA, Nexch)
for (i in 1:Nexch) {
  prior <- getPriorParameters("exnex", target_rates = pguess[i], tau_scale = tau.HN.scale[i])
  mu.mean[i] <- prior$exnex$mu_mean
  mu.prec[i] <- 1 / prior$exnex$mu_sd^2
}
prior <- getPriorParameters("exnex", target_rates = pguess[Nmix])
nex.mean <- prior$exnex$mu_j
nex.prec <- 1 / prior$exnex$tau_j^2

## start simulation
seed <- 887987
start_time <- Sys.time()
res.post <- post.infer.parallel(nclust = nclust, nperclust = nperclust, data.object = data.object, 
                                pnull = pnull, stopbounds =  stopbounds, beta.a0 = pnull, 
                                beta.b0 = 1-pnull, seed = seed, ModelFit = "EXNEX", 
                                weight = weight, nex.mean = nex.mean, nex.prec = nex.prec, 
                                mu.mean = mu.mean, mu.prec = mu.prec, tau.HN.scale = tau.HN.scale)
(Q <- get.Q.bwer(res.post, alpha = sig.level, digits = ndigits, Qclust = rep(1, B)))
res <- get.weighted.power(res.post, Q = Q, alpha = sig.level)
end_time <- Sys.time()
(time.EXNEX <- end_time - start_time)
## save results
results.EXNEX <- list(res.post, res)
Qmat <- array(NA, dim = dim(res.post$postprob))
for(i in 1:B) Qmat[,,i] <- res$Q[i]
(oc.EXNEX <- apply(res.post$postprob>Qmat, c(1,3), mean))
save.image("intermediate_results/AllMethods.RData")

####################### CBHM ########################
## prior values
phi = 0.5; # phi tuning for clustering
w = 2; # # w tuning for clustering
mu0 = 0; # hyperprior mean
sig02 = 1e6; # hyperprior variance                 
a0 = 1e-6; # hyperprior gamma(a0, b0) for 1/variance
b0 = 1e-6; # hyperprior gamma(a0, b0) for 1/variance
a1 = 0.1; # beta(a1, b1) prior for size one clusters.
b1 = 0.1; # beta(a1, b1) prior for size one clusters.
threshold_f <- (pnull+ptarget)/2
## start simulation
seed <- 887987
start_time <- Sys.time()
res.post <- post.infer.parallel(nclust = nclust, nperclust = nperclust, data.object = data.object, 
                                pnull = pnull, stopbounds =  stopbounds, beta.a0 = pnull, 
                                beta.b0 = 1-pnull, seed = seed, ModelFit = "CBHM", 
                                threshold_c = threshold_f, NDat = N[,ncol(N)],
                                phi = phi, w = w, mu0 = mu0, sig02 = sig02, a0 = a0, b0 = b0, 
                                a1 = a1, b1 = b1)
(Q <- get.Q.bwer(res.post, alpha = sig.level, digits = ndigits, Qclust = rep(1, B)))
res <- get.weighted.power(res.post, Q = Q, alpha = sig.level)
end_time <- Sys.time()
(time.CBHM <- end_time - start_time)
## save results
results.CBHM <- list(res.post, res)
Qmat <- array(NA, dim = dim(res.post$postprob))
for(i in 1:B) Qmat[,,i] <- res$Q[i]
(oc.CBHM <- apply(res.post$postprob>Qmat, c(1,3), mean))
save.image("intermediate_results/AllMethods.RData")

####################### MEM ########################
## prior values
shape1 <- pnull[1]
shape2 <- 1- pnull[1]
## start simulation
seed <- 887987
start_time <- Sys.time()
res.post <- post.infer.parallel(nclust = nclust, nperclust = nperclust, data.object = data.object, 
                                pnull = pnull, stopbounds =  stopbounds, beta.a0 = pnull, 
                                beta.b0 = 1-pnull, seed = seed, ModelFit = "MEM", 
                                shape1 = shape1, shape2 = shape2)
(Q <- get.Q.bwer(res.post, alpha = sig.level, digits = ndigits, Qclust = rep(1, B)))
res <- get.weighted.power(res.post, Q = Q, alpha = sig.level)
end_time <- Sys.time()
(time.MEM <- end_time - start_time)
## save results
results.MEM <- list(res.post, res)
Qmat <- array(NA, dim = dim(res.post$postprob))
for(i in 1:B) Qmat[,,i] <- res$Q[i]
(oc.MEM <- apply(res.post$postprob>Qmat, c(1,3), mean))
save.image("intermediate_results/AllMethods.RData")

