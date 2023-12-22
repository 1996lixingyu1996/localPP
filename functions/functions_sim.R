## Generate data
## for simulation replicates: trial=1, ..., ntrial;
##     indications: k=1, ..., B. 
##     scenarios: s=1, ..., nS;
## Returned object: array with dim = c(nS, ntrial, B, stage)
generate.data <- function(N, ORRs, ntrial = 10000, seed = 987897){
  # N: matrix with dim=(B, stage), where stage is the # of analyses (interim+final)
  # ORRs: a matrix with dim = (nS, B)
  set.seed(seed)
  if(is.vector(ORRs)){
    ORRs <- matrix(ORRs, nrow = 1)
  }
  if(is.vector(N)){
    N <- matrix(N)
  }
  stage <- dim(N)[2] # number of analyses (interim+final)
  B <- ncol(ORRs) # Number of Indications 
  nS <- nrow(ORRs) # number of scenarios
  res <- array(NA, dim = c(nS, ntrial, B, stage))  
  for (s in 1:nS){
    for (trial in 1:ntrial) {
      # Generate Observations for each Indication
      res[s, trial, , 1] <- rbinom(B, N[, 1], ORRs[s,])
      if(stage>1){
        for(j in 2:stage){
          res[s, trial, , j] <- res[s, trial, , j-1] + rbinom(B, N[, j]-N[, j-1], ORRs[s,])
        }
      }
    }
  }
  return(list(data = res, N = N, ORRs = ORRs))
}

## Generate posterior probabilities P(p_j > pnull_j) after all interim analysis
## and calculate rates for early stopping, number of patients and estimated ORR
## Returned object: array with dim = c(nS, ntrial, B)
post.infer <- function(object, pnull, stopbounds = NULL, clusterk = NULL, nperclust = NULL,
                       beta.a0 = pnull, beta.b0 = 1-pnull, seed = 987897, ModelFit, ...) {
  # object: returned from generate.data
  # pnull: B-vector of null response rates
  # stopbounds: B by (stage-1) matrix: stopping boundaries for each indication at each interim
  # beta.a0, beta.b0: B-vector of prior values for the beta prior for each indication
  ORRs <- object$ORRs
  N <- object$N
  if(is.null(clusterk) | is.null(nperclust)){
    dat <- object$data
  }else{
    dat <- object$data[, (clusterk-1)*nperclust+(1:nperclust), , ,drop=FALSE]
  }
  stage <- dim(N)[2] # number of analyses (interim+final)
  Nmax <- N[,stage]
  B <- ncol(ORRs) # Number of Indications 
  nS <- nrow(ORRs) # number of scenarios
  ntrial <- dim(dat)[2]
  Fit <- get(ModelFit)
  is.mc <- !(ModelFit %in% c("Independent", "localPP", "JSD", "localMEM"))
  is.MEM <- (ModelFit == "MEM")
  
  # Simulate Trials
  res.post <- array(NA, dim = c(nS, ntrial, B)) # posterior probabilities
  res.estop <- array(NA, dim = c(nS, ntrial, B)) # early stop status
  res.pts <- array(NA, dim = c(nS, ntrial, B)) # final enrolled number of patients 
  res.est <- array(NA, dim = c(nS, ntrial, B)) # estimate of ORR
  
  set.seed(seed)
  
  for (s in 1:nS){
    for (trial in 1:ntrial) {
      # read observations for each Indication
      yobs <- array( dat[s, trial, , ], dim = dim(dat)[3:4] )
      if(stage==1){
        y <- yobs[, stage]
        n <- N[,stage]
        if(is.mc){
          if(is.MEM){
            fit0 <- Fit(nDat = n, yDat = y, p0 = pnull[1], ...)
            pp <- fit0$post_prob
            phat <- fit0$phat
          }else{
            fit0 <- Fit(nDat = n, yDat = y, ...)
            pp <- colMeans(fit0>matrix(pnull, nrow(fit0), ncol(fit0), byrow = TRUE))
            phat <- colMeans(fit0)
          }
        }else{
          fit0 <- Fit(nDat = n, yDat = y, be.a0 = beta.a0, be.b0 = beta.b0, ...)
          pp <- pbeta(pnull, fit0$a.post, fit0$b.post, lower.tail = FALSE)
          phat <- fit0$a.post/(fit0$a.post+fit0$b.post)
        }
        res.post[s, trial, ] <- pp  
        res.estop[s, trial, ] <- rep(0, B)
        res.pts[s, trial, ] <- n
        res.est[s, trial, ] <- phat
      }else{
        Last_Stage = rep(1, B) # Keep Track of Last Stage
        y <- yobs[, 1]
        n <- N[, 1]
        stop.flag <- rep(FALSE, B)
        # interim analysis
        for (j in 1:(stage-1)){
          ind <- which(stop.flag==FALSE)
          if(length(ind)==0) break
          stop.flag[ind] <- (y[ind] <= stopbounds[ind,j])
          Last_Stage <- Last_Stage + !stop.flag
          y <- yobs[cbind(1:B, Last_Stage)] 
          n <- N[cbind(1:B, Last_Stage)]
        }
        res.post[s, trial, ] <- pbeta(pnull, beta.a0+y, beta.b0+n-y, lower.tail = FALSE)
        res.estop[s, trial, ] <-  stop.flag+0
        res.est[s, trial, ] <- ((beta.a0+y)/(beta.b0+n-y))
        res.pts[s, trial, ] <- n 
        
        # Final Analysis, where n contains the Final Sample Size
        ind.left <- which(stop.flag==FALSE)
        if(length(ind.left)>1){
          if(is.mc){
            if(is.MEM){
              fit0 <- Fit(nDat = n[ind.left], yDat = y[ind.left], p0 = pnull[1], ...)
              pp <- fit0$post_prob
              phat <- fit0$phat
            }else{
              fit0 <- Fit(nDat = n[ind.left], yDat = y[ind.left], ...)
              pp <- colMeans(fit0>matrix(pnull[ind.left], nrow(fit0), ncol(fit0), byrow = TRUE))
              phat <- colMeans(fit0)
            }
          }else{
            fit0 <- Fit(nDat = n[ind.left], yDat = y[ind.left], 
                        be.a0 = beta.a0[ind.left], be.b0 = beta.b0[ind.left], ...)
            pp <- pbeta(pnull[ind.left], fit0$a.post, fit0$b.post, lower.tail = FALSE)
            phat <- (fit0$a.post/(fit0$a.post+fit0$b.post))
          }
          res.post[s, trial, ind.left] <- pp
          res.est[s, trial, ind.left] <- phat
        }
      }
    }
  }
  return (list(earlystop = res.estop, postprob = res.post, npts = res.pts, est = res.est, 
               N = N, ORRs = ORRs, pnull = pnull, stopboundary = stopbounds))
}

## Get Q based on basket-wise error rate (bwer) control
get.Q.bwer <- function(object, alpha = 0.1, digits = 3, Qclust = NULL){
  #object: returned by post.infer()
  #Qclust: =NULL means all Qs are different;
  #        If there are B=5 baskets and Qclust=(1,1,2,2,2), it means Q for
  #        the first two baskets will be the same and another Q will be used for 
  #        baskets 3-5. 
  res.post <- object$postprob ### dim(nS, ntrial, B)
  pnull <- object$pnull
  Nmax <- object$N[,ncol(N)]
  nS <- dim(res.post)[1]
  ntrial <- dim(res.post)[2]
  B <- dim(res.post)[3]
  ORRs <- object$ORRs
  Q.final <- rep(NA, B)
  if(is.null(Qclust)) Qclust <- 1:B
  Q.unique <- unique(Qclust)
  nQ <- length(Q.unique)
  for(i in 1:nQ){
    i.ind <- which(Qclust==Q.unique[i])
    Q <- quantile(as.vector(res.post[1, , i.ind]), probs = 1-alpha, names = FALSE)
    Q <- ceiling(Q*10^digits)/(10^digits)
    Q.final[i.ind] <- ceiling(Q*10^digits)/(10^digits)
  }
  Q.final
}


## Get weighted type I error (WE) and power(WP) cross all scenarios.
## Including family wise (fwer) or trial wise (twer) or false discovery rate (fdr)
get.weighted.power <- function(object, Q, alpha = 0.1, s0 = 100, s1 = 0){
  # Setting s0=100 the weighted power reduces to type I error under global null
  # Setting s1=0 gives equal weight for calculating weighted power across scenarios.
  #object: returned by post.infer()
  res.post <- object$postprob ### dim(nS, ntrial, B)
  pnull <- object$pnull
  nS <- dim(res.post)[1]
  ntrial <- dim(res.post)[2]
  B <- dim(res.post)[3]
  ORRs <- object$ORRs
  pnull.mat <- matrix(pnull, nS, B, byrow = TRUE)
  ## weights
  H0status <- (ORRs<=pnull.mat) 
  nH0 <- rowSums(H0status)
  ind0 <- which(nH0!=0)
  wei0 <- rep(0, nS)
  wei0[ind0] <-  nH0[ind0]^s0/sum(nH0[ind0]^s0) 
  H1status <- (ORRs>pnull.mat) 
  nH1 <- rowSums(H1status)
  ind1 <- which(nH1!=0)
  wei1 <- rep(0, nS)
  wei1[ind1] <-  nH1[ind1]^s1/sum(nH1[ind1]^s1) 
  
  # Errors and powers
  Qarray <- array(NA, dim = c(nS, ntrial, B))  
  for(i in 1:B) Qarray[,,i] <- Q[i]
  res.rej <- (res.post>Qarray)
  fdrs <- rep(NA, nS)
  twerrors <- rep(NA, nS)
  fwerrors <- rep(NA, nS)
  cdrs <- rep(NA, nS)
  adrs <- rep(NA, nS)
  ccrs <- rep(NA, nS)
  bwer <- matrix(NA, nS, B)
  
  for(i in 1:nS){
    rate.rej <- matrix(res.rej[i,,], ntrial, B)
    x0 <- rate.rej[,ORRs[i,]<=pnull,drop=FALSE]
    x1 <- rate.rej[,ORRs[i,]>pnull,drop=FALSE]
    n.rej <- rowSums(rate.rej)
    n.fd <- ifelse(n.rej==0, 0, rowSums(x0)/n.rej) 
    fdrs[i] <- ifelse(length(x0)==0, NA, mean(n.fd))
    fwerrors[i] <- ifelse(length(x0)==0, NA, mean(apply(x0, 1, any)))
    twerrors[i] <- ifelse(length(x0)==0, NA, mean(colMeans(x0))) 
    #fwps[i] <- ifelse(length(x1)==0, NA, mean(apply(x1, 1, any)))
    #twps[i] <- ifelse(length(x1)==0, NA, mean(colMeans(x1))) 
    cdrs[i] <- ifelse(length(x1)==0, NA, mean(rowMeans(x1)))
    adrs[i] <- ifelse(length(x1)==0, NA, mean((rowSums(x1)-rowSums(x0))/ncol(x1)))
    ccrs[i] <- ifelse(length(x1)==0, NA, mean((rowSums(x1)+(ncol(x0)-rowSums(x0)))/B)) 
    bwer[i,] <- colMeans(rate.rej)
    bwer[i,ORRs[i,]>pnull] <- NA
  }
  w.fdr <- weighted.mean(fdrs, w = wei0)
  w.twe <- weighted.mean(twerrors, w = wei0)
  w.fwe <- weighted.mean(fwerrors, w = wei0)
  w.cdr <- weighted.mean(cdrs, w = wei1)
  w.adr <- weighted.mean(adrs, w = wei1)
  w.ccr <- weighted.mean(ccrs, w = wei1)
  
  return(list(Q=Q, error.fdr = w.fdr, error.tw=w.twe, error.fw=w.fwe, 
              power.cdr=w.cdr, power.adr = w.adr, power.ccr = w.ccr,
              ind.error.fdr = fdrs, ind.error.tw = twerrors, ind.error.fw = fwerrors, 
              ind.power.cdr=cdrs, ind.power.adr=adrs, ind.power.ccr=ccrs, bwer = bwer,
              w0 = wei0, w1 = wei1))
}

