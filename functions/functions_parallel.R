## Generate posterior probabilities parallel version
post.infer.parallel <- function(nclust, nperclust, data.object, pnull, stopbounds = NULL, 
                                beta.a0 = pnull, beta.b0 = 1-pnull, seed = 38098, ModelFit, ...) {
  if(nclust > detectCores()) stop("nclust needs to be smaller than detectCores()")
  
  ## Create Clusters
  cl <- makeCluster(nclust, type = "PSOCK")
  registerDoParallel(cl)
  
  ## Export Functions to the Cluster
  clusterCall(cl, function() {source("../functions/functions.R"); source("../functions/functions_sim.R")})
  clusterExport(cl, list("nperclust", "data.object", "pnull", "stopbounds", "beta.a0", "beta.b0", "seed"))
  
  ## Parallel Computing
  results <- parLapply(cl, (1:nclust), fun = post.infer, object = data.object, 
                       pnull = pnull, stopbounds = stopbounds, nperclust = nperclust, 
                       beta.a0 = beta.a0, beta.b0 = beta.b0, seed = seed,
                       ModelFit = ModelFit, ...)
  
  ## Close Cluster
  stopCluster(cl)
  
  ## Combine the Results
  resultsAll <- results[[1]]
  for(i in 2:nclust){
    resultsAll$earlystop <- abind::abind(resultsAll$earlystop, results[[i]]$earlystop, along=2)
    resultsAll$postprob <- abind::abind(resultsAll$postprob, results[[i]]$postprob, along=2)
    resultsAll$npts <- abind::abind(resultsAll$npts, results[[i]]$npts, along=2)
    resultsAll$est <- abind::abind(resultsAll$est, results[[i]]$est, along=2)
  }
  return(resultsAll)
}

