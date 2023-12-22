
### Upload Required Packages
library(parallel)
library(doParallel)
library(R2jags)
library(foreach)
library(bhmbasket) ## required for EXNEX
library(basket) ## required for MEM (Kane et al., 2020)
source("../functions/utils.R") ## required for BCHM (Chen and Lee, 2020)
source("../functions/cluster.R") ## required for local MEM (Liu et al. 2022)


### logit function
logit<-function(x)
{
  log(x/(1-x))
}

### Independent beta binomial model
Independent <- function(nDat, yDat, be.a0 = NULL, be.b0 = NULL){
  B <- length(nDat)
  if(is.null(be.a0) | is.null(be.b0)) {
    be.a0 = rep(.15, B)
    be.b0 = rep(.85, B)
  }
  return(list(a.post = be.a0 + yDat, b.post = be.b0 + nDat - yDat))
}

## MEM (Kane et al., 2020)
MEM <- function(nDat, yDat, p0 = 0.15, shape1 = 0.15, shape2 = 0.85){
  n <- nDat
  y <- yDat
  B <- length(n)
  res <- mem_exact(responses = y, size = n, name = letters[1:B], p0 = p0, 
                   shape1 = shape1,  shape2 = shape2)
  return(list(phat = as.vector(res$basket$mean_est), 
              post_prob = as.vector(res$basket$post_prob)))
}

## local Power Prior
localPP <- function(nDat, yDat, be.a0 = NULL, be.b0 = NULL, a = 0.4,
                    delta = 0.15, symmetry = FALSE){
  n <- nDat
  y <- yDat
  B <- length(n)
  if(is.null(be.a0) | is.null(be.b0)) {
    be.a0 = rep(.15, B)
    be.b0 = rep(.85, B)
  }
  if(B==1){
    return(list(a.post = be.a0 + y, b.post = be.b0 + n - y, sm = matrix(1)))
  }else{
    # similarity matrix
    sm <- matrix(1, B, B)
    fn <- function(theta, yi, ni, y0, n0, a0, b0) {
      (lbeta(a0+yi+theta*y0, b0+(ni-yi)+theta*(n0-y0)) 
       -lbeta(a0+theta*y0, b0+theta*(n0-y0)))
    }
    for(i in 1:B){
      for(j in (1:B)[-i]){
        res <- optimize(f = fn, interval = c(0, 1),  
                        yi = y[i], ni = n[i], y0 = y[j], n0 = n[j], 
                        a0 = be.a0[i], b0 = be.b0[i], maximum = TRUE)
        sm[i,j] <- res$maximum*a*(abs(y[i]/n[i]-y[j]/n[j])<delta)
      }
    }
    if (symmetry) sm <- (t(sm) + sm)/2
    ab.post <- matrix(NA, B, 2)
    for(i in 1:B){
      ab.post[i,1] <- be.a0[i]+sum(sm[i,]*y)
      ab.post[i,2] <- be.b0[i]+sum(sm[i,]*(n-y))
    }
    return(list(a.post = ab.post[,1], b.post = ab.post[,2], 
                sm = sm))
  }
}

## JSD (Fujikawa et al., 2020)
JSD <- function(nDat, yDat, be.a0 = NULL, be.b0 = NULL, epsilon = 2, tau = 0.3){
  n <- nDat
  y <- yDat
  B <- length(n)
  if(is.null(be.a0) | is.null(be.b0)) {
    be.a0 = rep(.15, B)
    be.b0 = rep(.85, B)
  }
  if(B==1){
    return(list(a.post = be.a0 + y, b.post = be.b0 + n - y))
  }else{
    # similarity matrix
    sm <- matrix(1, B, B)
    for(i in 1:(B-1)){
      fi <- function(x) dbeta(x, be.a0[i]+y[i], be.b0[i] + n[i] - y[i])
      for(j in (i+1):B){
        fj <- function(x) dbeta(x, be.a0[j]+y[j], be.b0[j] + n[j] - y[j])
        tmp1 <- integrate(function(x) fi(x)*log(2*fi(x)/(fi(x)+fj(x))), 
                          lower=0.0001, upper=0.9999 )$value
        tmp2 <- integrate(function(x) fj(x)*log(2*fj(x)/(fi(x)+fj(x))), 
                          lower=0.0001, upper=0.9999 )$value
        smij <- 1-(tmp1+tmp2)/2
        #sm[i,j] <- sm[j,i] <- smij*a*(abs(yo[i]/no[i]-yo[j]/no[j])<=delta)
        #pij <- prop.test(yo[c(i,j)], no[c(i,j)])$p.value
        #sm[i,j] <- sm[j,i] <- smij*a*(pij>=delta)
        sm[i,j] <- sm[j,i] <- (smij^epsilon)*((smij^epsilon)>tau)
      }
    }
    ab.post <- matrix(NA, B, 2)
    for(i in 1:B){
      ab.post[i,1] <- be.a0[i]+sum(sm[i,]*y)
      ab.post[i,2] <- be.b0[i]+sum(sm[i,]*(n-y))
    }
    return(list(a.post = ab.post[,1], b.post = ab.post[,2]))
  }
}

## local MEM (Liu et al., 2022)
localMEM <- function(nDat, yDat, be.a0 = NULL, be.b0 = NULL, 
                     delta = 2, # Delta Specification for the Prior
                     max_C = length(nDat) # Max # of Clusters
                     ){
  n <- nDat
  y <- yDat
  B <- length(n)
  if(is.null(be.a0) | is.null(be.b0)) {
    be.a0 = rep(.15, B)
    be.b0 = rep(.85, B)
  }
  if(B==1){
    return(list(a.post = be.a0 + y, b.post = be.b0 + n - y))
  }else{
    ## Partition
    part <- get.part(R = B, max_cl = max_C)
    K <- nrow(part)
    # Number of Blocks/Unique Response Rate in Each Partition
    n_bk <- apply(part, 1, function(x){
      length(unique(x))
    })
    prior_part <- n_bk^delta/sum(n_bk^delta)
    post_part <- update.part(x = y, n = n, prior_part = prior_part, 
                             part = part, a0 = 0.5, b0 = 0.5)
    fpart <- unlist( part[which.max(post_part), ])
    fmat <- NULL
    for(r in 1:B){
      fmat <- rbind(fmat, ifelse(fpart[r]==fpart, 1, 0))
    }
    diag(fmat) <- 0
    pOmega <- max(post_part)
    ab.post <- matrix(NA, B, 2)
    for(i in 1:B){
      ab.post[i,1] <- be.a0[i]+y[i]+sum(fmat[i,]*y)*pOmega
      ab.post[i,2] <- be.b0[i]+n[i]-y[i]+sum(fmat[i,]*(n-y))*pOmega
    }
    return(list(a.post = ab.post[,1], b.post = ab.post[,2]))
  }
}


# BHM uniform (Cunanan et al., 2019)
BHMunif <- function(nDat,
                yDat,
                mu0 = 0, # hyperprior for mu~N(mu0, sig02)
                sig02 = 100, #hyperprior for mu~N(mu0, sig02)                 
                u0 = 100, # hyperprior for sig~ U(0, u0) 
                be.a0 = NULL, # Beta(a0, b0) on ORR for each basket for clustering
                be.b0 = NULL, # Beta(a0, b0) on ORR for each basket for clustering
                n.chains = 2,
                n.iter = 12000, 
                n.burnin = 2000,
                n.thin = 5)
{
  n <- nDat
  y <- yDat
  B <- length(n)
  nsave <- (n.iter-n.burnin)/n.thin*n.chains
  
  if(is.null(be.a0) | is.null(be.b0)) {
    be.a0 = rep(.15, B)
    be.b0 = rep(.85, B)
  }
  
  psamples <- matrix(NA, nsave, B)
  if(B==1){
    psamples[,1] <- rbeta(nsave, be.a0+y, be.b0+n-y)
  }else{
    jagsmodel <- function()
    {
      for (j in 1:B){
        y[j] ~ dbinom(p[j], n[j]) # Likelihood
        logit(p[j]) <- theta[j]
        theta[j] ~ dnorm(mu, tau2) 
      }
      mu ~ dnorm(mu0, 1/sig02)         # Prior on mu
      sig ~ dunif(0, u0)    # Prior on sig
      tau2 <- 1/sig^2
    }
    # Keep track of p
    jags.params<-c("p")
    # Initial Value for Sampling Scheme
    jags.inits <-function(){
      list("mu" = 0, "sig" = 1)
    }
    jags.data <- list("B", "y", "n", "mu0", 
                      "sig02", "u0")
    jagsfit <- jags(data = jags.data, inits = jags.inits, parameters.to.save = jags.params,
                    n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin,
                    model.file = jagsmodel, quiet = TRUE, progress.bar = "none")
    psamples <- jagsfit$BUGSoutput$sims.matrix[,-1, drop = FALSE]
  }
  return(psamples)
}

# BHM gamma prior (Berry et al., 2013)
BHM <- function(nDat,
                yDat,
                mu0 = -1.735, # Hyperprior Mean
                sig02 = 10, # Hyperprior Variance                 
                a0 = 2, # Hyperprior Gamma(a0, b0) for 1/Variance
                b0 = 20, # Hyperprior Gamma(a0, b0) for 1/Variance
                n.chains = 2,
                n.iter = 12000, 
                n.burnin = 2000,
                n.thin = 5)
{
  n <- nDat
  y <- yDat
  numGroups <- length(n)
  
  MeanMu <- mu0       
  PercMu <- 1/sig02  
  TauAlpha <- a0
  TauBeta <- b0
  
  jagsmodel <- function()
  {
    mu ~ dnorm(MeanMu, PercMu)         # Prior on Mu
    tau ~ dgamma(TauAlpha, TauBeta)    # Prior on Tau
    
    for (j in 1:numGroups){
      y[j] ~ dbinom(p[j], n[j]) # Likelihood
      theta[j] ~ dnorm(mu, tau) 
      p[j] <- 1/(1 + exp(-theta[j]))
    }
    
  }
  # Keep track of p
  jags.params<-c("p")
  # Initial Value for Sampling Scheme
  mu.init <- (sum(y)/sum(n) + 0.01)
  jags.inits <-function(){
    list("mu" = ifelse(mu.init < 0.99, logit(mu.init), 4.59512), "tau" = 1)
  }
  jags.data <- list("numGroups", "y", "n", "MeanMu", 
                    "PercMu", "TauAlpha", "TauBeta")
  jagsfit <- jags(data = jags.data, inits = jags.inits, parameters.to.save = jags.params,
                  n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin,
                  model.file = jagsmodel, quiet = TRUE, progress.bar = "none")
  return(jagsfit$BUGSoutput$sims.matrix[,-1, drop = FALSE])
}


## CBHM (Jiang et al. 2021) 
CBHM <- function(nDat,
                 yDat,
                 threshold_c, # (p.null + p.target)/2
                 NDat, # a vector of Total Sample Size for Each Indication >= nDat
                 phi = 0.5, # phi Tuning for Clustering
                 w = 2, # Tuning for Clustering
                 mu0 = 0, # Hyperprior Mean
                 sig02 = 1e6, # Hyperprior Variance                 
                 a0 = 1e-6, # Hyperprior Gamma(a0, b0) for 1/Variance
                 b0 = 1e-6, # Hyperprior Gamma(a0, b0) for 1/Variance
                 a1 = 0.1, # Beta(a1, b1) Prior for Size One Clusters.
                 b1 = 0.1, # Beta(a1, b1) Prior for Size One Clusters.
                 n.chains = 2,
                 n.iter = 12000, 
                 n.burnin = 2000,
                 n.thin = 5)
{
  n <- nDat
  y <- yDat
  N <- NDat
  numGroups <- length(n)
  nsave <- (n.iter-n.burnin)/n.thin*n.chains
  
  psamples <- matrix(NA, nsave, numGroups)
  if(numGroups==1){
    psamples[,1] <- rbeta(nsave, a1 + y, b1 + n - y)
  }else{
    # Clustering Indications
    eff_ind <- NULL
    for(g in 1:numGroups){
      if(pbeta(threshold_c[g], a1 + y[g], b1 + n[g] - y[g], lower.tail = FALSE) > phi*(n[g]/N[g])^w) { # Computes Which Clusters are "Effective"
        eff_ind <- c(eff_ind, g)
      }
    }
    all_ind <- (1:numGroups)
    non_ind <- setdiff(all_ind, eff_ind)
    
    if (length(eff_ind)==1){
      psamples[,eff_ind] <- rbeta(nsave, a1 + y[eff_ind], b1 + n[eff_ind] - y[eff_ind])
    } else if(length(eff_ind)>=2){
      psamples[,eff_ind] <- BHM(nDat = n[eff_ind], yDat = y[eff_ind], mu0 = mu0, sig02 = sig02, a0 = a0,  
                                b0 = b0, n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin,
                                n.thin = n.thin)
    }
    if(length(non_ind)==1){
      psamples[,non_ind] <- rbeta(nsave, a1 + y[non_ind], b1 + n[non_ind] - y[non_ind])
    }else if(length(non_ind)>=2){
      psamples[,non_ind] <- BHM(nDat = n[non_ind], yDat = y[non_ind], mu0 = mu0, sig02 = sig02, a0 = a0, b0 = b0, n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin,
                                n.thin = n.thin)
    }
  }
  return(psamples)
}

## BCHM (Chen and Lee 2020)
### This function generates samples from BCHM.
BCHM <- function(nDat, yDat,
                 mu = 0.2,
                 sigma02 = 10, # Variance of Prior
                 sigmaD2 = 0.001, # Variance of Each Data Point 
                 alpha = 1e-40, # Similarity Constant
                 d0 = 0, # Lower Bound on Similarity Matrix
                 alpha1 = 50, # for Gamma(alpha1, beta1) in BHM
                 beta1 = 10, # for Gamma(alpha1, beta1) in BHM
                 mu0 = -1.735, # for the Hyperprior Mean Mean in BHM
                 tau2 = 0.1, # for the Hyperprior Mean in BHM Variance/Precision
                 n.chains = 2,
                 n.iter = 12000,
                 n.burnin = 2000,
                 n.thin = 5) {
  numArm <- length(nDat)
  nsave <- (n.iter-n.burnin)/n.thin*n.chains
  if (numArm != length(yDat)) {
    stop("Numbers of subgroups in nDat and yDat are not equal.")
  }
  if (numArm > 20) {
    stop("Numbers of subgroups is more than 20.")
  }
  if (numArm==1){
    allPost <- BHM(nDat = nDat, yDat = yDat, mu0 = mu0, sig02 = 1/tau2, a0 = alpha1,  
                   b0 = beta1, n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin,
                   n.thin = n.thin)
  } else{
    weight <- nDat
    alphaP <- alpha
    alpha <- alpha1
    beta <- beta1
    posi <- yDat
    priorMean <- mu
    priorVar <- sigma02
    res <- posi/weight
    x <- as.matrix(res)
    estVarGroup <- sigmaD2
    result <- gibbsSampler(x, alphaP, priorMean, priorVar, estVarGroup, 
                           weight, 1000, 2000)
    tables <- result$tables
    sm <- matrix(0, numArm, numArm)
    tSize <- dim(tables)[1]
    for (i in 1:tSize) {
      rr <- tables[i, ]
      for (j in 1:numArm) {
        for (k in 1:numArm) {
          if (rr[j] == rr[k]) {
            sm[j, k] <- sm[j, k] + 1
          }
        }
      }
    }
    sm <- sm/tSize
    smMinR <- d0
    sm[sm < smMinR] <- smMinR
    SMatrix <- sm
    jagsmodel <- function()
    {
      for (i in 1:numGroups)
      {
        y[i] ~ dbin(p[i],n[i]) 
        logit(p[i]) <- eta[i] 
        rr[i]<-tau1 * m[i]
        eta[i] ~ dnorm(mu,rr[i])
      }
      # Priors
      mu ~ dnorm(mu0, tau2)
      tau1 ~ dgamma(alpha, beta)
    }
    
    if(all(sm > 0.99)){
      allPost <- BHM(nDat = weight, yDat = posi, mu0 = mu0, sig02 = 1/tau2, a0 = alpha1,  
                     b0 = beta1, n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin,
                     n.thin = n.thin)
    } else {
      for (i in 1:numArm) {
        smP <- sm[i, ]
        ind <- which(smP>0.01)
        n.ind <- length(ind)
        mydata <- list(y = posi[ind], n = weight[ind], m = smP[ind], 
                       numGroups = n.ind, mu0 = mu0, 
                       tau2 = tau2, alpha = alpha, beta = beta)
        # Parameters of Interest
        jags.params<-c("p")
        # Initial Value 
        mu.init <- sum(posi[ind])/sum(weight[ind]) + 0.01
        jags.inits<-function(){
          list("mu" = ifelse(mu.init<0.99, logit(mu.init), 4.59512), "tau1" = 1)
        }
        jagsfit <- jags(data = mydata, inits = jags.inits, parameters.to.save = jags.params, n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, model.file = jagsmodel, quiet = TRUE, progress.bar = "none")
        samples <- jagsfit$BUGSoutput$sims.matrix[,-1, drop = FALSE]
        sampledP <- samples[, which(ind==i)]
        
        if (i == 1) {
          allPost <- sampledP
        }
        else {
          allPost <- cbind(allPost, sampledP)
        }
      }
    }
  }
  allPost
}

## EXNEX 2015 Neuschwander Paper Recommendations
### This function generates samples from EXNEX.
EXNEX <- function(nDat, yDat, weight = c(0.25, 0.25, 0.5), 
                  nex.mean = -1.73, # NeX.Mean
                  nex.prec = 0.09, # Nex.Prec
                  mu.mean = c(-2.20, -2.20), # Mu.Mean in EX
                  mu.prec = c(0.099, 0.099), # Mu.Prec in EX
                  tau.HN.scale = c(1,1), # from the Paper - Prior Tau Precision Power Parameter
                  n.chains = 2,
                  n.iter = 12000,
                  n.burnin = 2000,
                  n.thin = 5) {
  Nexch = length(mu.mean)
  Nmix = Nexch + 1
  
  EXNEX_Model <- function(){
    
    # Prior Distribution for EX-parameters
    for (jj in 1:Nexch){
      mu[jj] ~ dnorm(mu.mean[jj], mu.prec[jj])
      prior.tau.prec[jj] <- pow(tau.HN.scale[jj],-2)
      tau[jj] ~ dnorm(0, prior.tau.prec[jj]); T(0.001, )
      prec.tau[jj] <- pow(tau[jj],-2)
    }
    
    # Log-Odds Parameters under NEX
    for (j in 1:Nstrata){
      LogOdds[Nmix,j] ~ dnorm(nex.mean, nex.prec)
    }
    
    # Log-Odds Parameters under EX
    for (jj in 1:Nexch){
      for (j in 1:Nstrata){
        re[jj,j] ~ dnorm(0, prec.tau[jj])
        LogOdds[jj,j] <- mu[jj]+re[jj,j]
      }
    }
    
    for (j in 1:Nstrata){
      exch.index[j] ~ dcat(pMix[1:Nmix])
    }
    
    for (j in 1:Nstrata){
      theta[j] <- LogOdds[exch.index[j],j]
    }
    
    for (i in 1:Nstrata){                
      logit(p[i]) <- theta[i]
      y[i] ~ dbin(p[i],n[i])
    }
  }
  
  # Initial Values
  jags.params <- c("p")
  
  jags.inits <- function(){
    list(mu = mu.mean)
  }
  
  Nstrata <- length(nDat)
  n <- nDat
  y <- yDat
  pMix <- weight
  
  jags.data <- list("y","n","Nstrata","Nexch","Nmix","nex.mean","nex.prec", "pMix", "mu.mean", "mu.prec", "tau.HN.scale")
  
  jagsfit <- jags(data = jags.data, inits = jags.inits, parameters.to.save = jags.params,
                  n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, model.file = EXNEX_Model, quiet = TRUE, progress.bar = "none")
  
  return (jagsfit$BUGSoutput$sims.matrix[,-1, drop = FALSE])
}

