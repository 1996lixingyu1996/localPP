####################################################################################
## Create the table for local PP with different hyperparameter a
####################################################################################
load("intermediate_results/localPP_select_a.RData")
results <- matrix(NA, length(aseq), 5)
rownames(results) <- paste("a=", aseq, sep="")
colnames(results) <- c("BWER_S1", "BWER-avg", "BWER-max", "TPR-avg", "ccr-avg")
for(i in 1:length(aseq)){
  res <- list.metrics0[[i]]
  results[i, ]<- sprintf(c(res$error.tw, mean(res$bwer, na.rm = TRUE), max(res$bwer, na.rm = TRUE),
                           res$power.cdr, res$power.ccr), fmt = '%#.3f')
}
results

####################################################################################
## Create the table for local PP based on pairwise similarity with different hyperparameter a
####################################################################################
load("intermediate_results/localPP_select_a_pairwise.RData")
results <- matrix(NA, length(aseq), 5)
rownames(results) <- paste("a=", aseq, sep="")
colnames(results) <- c("BWER_S1", "BWER-avg", "BWER-max", "TPR-avg", "ccr-avg")
for(i in 1:length(aseq)){
  res <- list.metrics0[[i]]
  results[i, ]<- sprintf(c(res$error.tw, mean(res$bwer, na.rm = TRUE), max(res$bwer, na.rm = TRUE),
                           res$power.cdr, res$power.ccr), fmt = '%#.3f')
}
results

####################################################################################
## Create the table for JSD with different hyperparameter tau
####################################################################################
load("intermediate_results/JSD_select_tau.RData")
results <- matrix(NA, length(aseq), 5)
rownames(results) <- paste("tau=", aseq, sep="")
colnames(results) <- c("BWER_S1", "BWER-avg", "BWER-max", "TPR-avg", "ccr-avg")
for(i in 1:length(aseq)){
  res <- list.metrics0[[i]]
  results[i, ]<- sprintf(c(res$error.tw, mean(res$bwer, na.rm = TRUE), max(res$bwer, na.rm = TRUE),
                           res$power.cdr, res$power.ccr), fmt = '%#.3f')
}
results

####################################################################################
## Create the rest of the tables
####################################################################################
load("intermediate_results/AllMethods.RData")
nmethods <- 8
method.names <- c("Independent", "local PP", "JSD", "EXNEX", "BHM", "BCHM", "local MEM", "MEM")

##### Table 3.2 The efficacy cutoff Q_is for all considered methods.
ii <- 1; 
res.post.all <- list(results.Independent[[ii]], results.localPP[[ii]],
                     results.JSD[[ii]], results.EXNEX[[ii]], results.BHMunif[[ii]], results.BCHM[[ii]], 
                     results.localMEM[[ii]], results.MEM[[ii]])
Q.by.method <- matrix(NA, nmethods, B+1)
rownames(Q.by.method) <- method.names
colnames(Q.by.method) <- c("Basket1", "Basket2", "Basket3", "Basket4", "Basket5", "Overall")
for(i in 1:nmethods){
  res.post.i <- res.post.all[[i]]
  Q.indi <- get.Q.bwer(res.post.i, alpha = sig.level, digits = ndigits, Qclust = 1:B)
  Q.overall <- get.Q.bwer(res.post.i, alpha = sig.level, digits = ndigits, Qclust = rep(1,B))
  Q.by.method[i,] <- sprintf(c(Q.indi, Q.overall[1]), fmt = '%#.3f')
}
Q.by.method

##### Table 3.3 Overall performance for different methods
results <- matrix(NA, nmethods, 6)
rownames(results) <- method.names
colnames(results) <- c("BWER_S1", "BWER-avg", "BWER-max", "TPR-avg", "CCR-avg", "Time (hours)")
res.oc.all <- list(results.Independent[[2]], results.localPP[[2]],
                    results.JSD[[2]], results.EXNEX[[2]], results.BHMunif[[2]], results.BCHM[[2]], 
                    results.localMEM[[2]], results.MEM[[2]])
time.all <- c(time.Independent, time.localPP, time.JSD, time.EXNEX, time.BHMunif, time.BCHM, 
                  time.localMEM, time.MEM)/3600
for(i in 1:nmethods){
  res <- res.oc.all[[i]]
  results[i, ]<- sprintf(c(res$error.tw, mean(res$bwer, na.rm = TRUE), max(res$bwer, na.rm = TRUE),
                           res$power.cdr, res$power.ccr, time.all[i]), 
                         fmt = '%#.3f')
}
results

## Table 3.4  Results all considered methods by scenarios
oc.allmethods <- list(oc.Independent, oc.localPP, oc.JSD, oc.EXNEX, oc.BHMunif, 
                      oc.BCHM, oc.localMEM, oc.MEM)
oc.by.sc <- matrix(NA, nmethods, 9)

rownames(oc.by.sc) <- method.names
colnames(oc.by.sc) <- c("Basket1", "Basket2", "Basket3", "Basket4", "Basket5", "FPR", "FDR", "TPR", "CCR")
res.by.sc <- list(scenario1=NULL, scenario2=NULL, scenario3=NULL,
                  scenario4=NULL, scenario5=NULL, scenario6=NULL)
for(k in 1:nrow(scenarios)){
  for(i in 1:nmethods){
    res <- res.oc.all[[i]]
    oc.i <- oc.allmethods[[i]]
    oc.by.sc[i,1:5]<- sprintf(oc.i[k,], fmt = '%#.3f')
    oc.by.sc[i,6:9]<- sprintf(c(res$ind.error.tw[k], res$ind.error.fdr[k], 
                                res$ind.power.cdr[k], res$ind.power.ccr[k]), fmt = '%#.3f')
  }
  res.by.sc[[k]] <- oc.by.sc
}
res.by.sc