## These functions are required for fitting BCHM (Chen and Lee, 2020)
library(plyr)
#library(clValid)
library(cluster)
library(coda)
library(rjags)

#' @importFrom cluster silhouette
#' @importFrom stats dist
#' 

modelStr<-function()
{
  mod1 <- "model
{
	for (i in 1:numGroups)
	{
		y[i] ~ dbin(p[i],n[i]);
		logit(p[i]) <- eta[i];
                rr[i]<-tau1 * m[i]
		eta[i] ~ dnorm(mu,rr[i])

		pg[i] <- step(p[i] - targetResp)
		# Probability that the response rate for each group is greater than targetResp	
	}
	#Priors
	mu ~ dnorm(mu0, tau2)
	tau1 ~ dgamma(alpha, beta)
}"
  return(mod1)
}


sigmoid<-function(x)
{
  1.0/(1+exp(-x))
}

logit<-function(x)
{
  log(x/(1-x))
}

varBeta<-function(alpha, beta)
{
  t<-alpha * beta / (alpha + beta) / (alpha + beta) /(alpha + beta + 1)
}

boa.hpd <- function(x, alpha) {
  n <- length(x)
  m <- max(1, ceiling(alpha * n))
  y <- sort(x)
  a <- y[seq_len(m)]
  b <- y[(n - m + 1):n]
  i <- order(b - a)[1]
  list(lower=a[i], upper=b[i])
}

getVar<-function(inSD)
{
  meanResp <- 0.2
  armSi <- 25
  al<- armSi * meanResp + 1
  be<- armSi - armSi * meanResp + 1
  sd2 <- 2 * sqrt(varBeta(al,be))    
  return(1+exp(20*(inSD-sd2)))
}


plotResult<-function(x, tab)
{
  plot(1:dim(x)[1], x, col=tab)
}


# ScaleFac
scaleFac <- function(allData) {
  
  renorm <- function(x,y)
  {
    if ((y == -Inf) && (x == -Inf))
    {
      return(-Inf)
    }
    
    if (y < x){
      return(x+log(1 + exp(y-x)))
    }  else {
      return(y+log(1 + exp(x-y)))
    }
  }
  
  if (length(allData) == 1){
    return(allData)
  }
  
  r <- allData[1];
  
  for (i in 2:length(allData)){
    r <- renorm(r, allData[i])
  }
  return(r)
}

## The distL and the sample functions are based on an example given by David Blei 

distL <- function(sumV, sumsq, num, pMean, pVar, dVar)
{
  varV <- 1 / (1 / pVar + num / dVar)
  meanV <- ((pMean / pVar) + (sumV / dVar)) * varV
  
  te <- - sumsq / dVar
  te <- te + (meanV * meanV) / varV + log(varV)
  te <- te - (pMean * pMean) / pVar - log(pVar)
  te <- (length(sumV) / 2) * te
  list(mean=meanV, var=varV, lhood=te)
}


sss<-function(ss, point, point.sq, weight, prior.mean, prior.var, data.var) {
  distL(ss$sum + point * weight,
            ss$sumsq + point.sq * weight,
            ss$count + 1 * weight,
            prior.mean, prior.var, data.var)$lhood - ss$lhood
}


llNew<-function(t, point,  point.sq, weight, prior.mean, prior.var, data.var)
{
  res<-c()
  for (i in 1:length(t))
  {
    res<-c(res, sss(t[[i]], point, point.sq, weight, prior.mean, prior.var, data.var))
  }
  res
}


sampleGroups <- function(point, allGroups, counts, alpha,
                         prior.mean, prior.var, data.var, weight)
{

  n <- sum(counts)
  p <- length(point)
  point.sq <- point %*% point
  

  log.denom <- log(n + alpha)

  log.prior <- c(log(counts) - log.denom, log(alpha) - log.denom)

  ll<-llNew(allGroups, point, point.sq, weight, prior.mean, prior.var, data.var)

  newL <- distL(point * weight, point.sq * weight, 1 * weight,
                    prior.mean, prior.var, data.var)$lhood
  
  prob <- log.prior + c(ll, newL)
  
  prob <- prob - scaleFac(prob)
  prob <- exp(prob)
  sample(1:length(prob), 1, prob=prob)
}


gibbsSampler <- function(data, alpha, prior.mean, prior.var, data.var, weight, burnIn=100, MCIter=200)
{
  
  #browser()
  n <- dim(data)[1]
  
  tables <- rep(0,n)
  groups <- list()
  counts <- c()
  allTab<-matrix(0, MCIter, n)
  
  for (iter in 1:(burnIn+MCIter))
  {
    #cat("iter = ", iter, "\n")
    for (i in 1:n)
    {
      # get table assignment and data point
      z <- tables[i]
      x <- data[i,]
      w <- weight[i]
      # remove this data point from its table
      if (z > 0)
      {
        groups[[z]]$sum <- groups[[z]]$sum - x * w
        groups[[z]]$sumsq <- groups[[z]]$sumsq - (x %*% x) * w
        groups[[z]]$count <- groups[[z]]$count - w
        groups[[z]]$lhood <- distL(groups[[z]]$sum, groups[[z]]$sumsq,
                                       groups[[z]]$count, prior.mean,
                                       prior.var, data.var)$lhood
        counts[z] <- counts[z] - w
        
        # Remove the group if number of points is 0
        if (counts[z] == 0)
        {
          # remove the table from the restaurant configuration
          counts <- counts[-z]
          groups <- groups[-z]
          # renumber the tables "beyond" this table
          tables[tables > z] <- tables[tables > z] - 1
        }
      }
      
      # Chinese restaurant process
      if (is.null(counts))
        z <- 1
      else
        z <- sampleGroups(x, groups, counts, alpha,
                          prior.mean, prior.var, data.var, w)
      
      tables[i] <- z
      
      if (z > length(counts))
      {
        
        stopifnot(z==length(counts) + 1)
        counts <- c(counts, 0)
        groups[[z]] <- list(sum=0, sumsq=0, lhood=0, count = 0)
      }
      
      counts[z] <- counts[z] + w
      groups[[z]]$sum <- groups[[z]]$sum + x * w
      groups[[z]]$sumsq <- groups[[z]]$sumsq + x  %*%  x * w
      groups[[z]]$count <- groups[[z]]$count + w
      groups[[z]]$lhood <- distL(groups[[z]]$sum, groups[[z]]$sumsq,
                                     groups[[z]]$count, prior.mean,
                                     prior.var, data.var)$lhood
    }

    if (iter > burnIn)
    {
      allTab[iter-burnIn, ]<-tables
    }
  }
  list(counts=counts, groups=groups, tables=allTab)
}


initData<-function(count, numGroup,sdGroup)
{
  center<-runif(numGroup, 0.1, 0.5)
  sdV<-sdGroup
  #center<-runif(numGroup, -10, 10)
  #var<-0.5
  dat<-c()
  for (i in 1: numGroup)
  {
    dat<-c(dat, rnorm(count, mean=center[i], sd=sdV))  
  }
  dat<-sample(dat)
  x<-matrix(0, count*numGroup, 1)
  x[,1]<-dat
  x
}

initDataWithWeight<-function(count, numGroup, sdGroup, distStart=0.1, distEnd=0.7, minSize=10, maxSize=25)
{
  center<-runif(numGroup, distStart, distEnd)
  sdV<-sdGroup
  dat<-c()
  weight<-c()
  groupID<-c()
  for (i in 1: numGroup)
  {
    dat<-c(dat, rnorm(count, mean=center[i], sd=sdV))  
    groupID<-c(groupID, rep(i,count))
  }
  #dat<-sample(dat)
  l<-length(dat)
  weight<-round(runif(l,minSize, maxSize))
  pos<-weight*dat
  pos<-round(pos)
  dat<-pos/weight
  
  x<-matrix(0, count*numGroup, 4)
  x[,1]<-dat
  x[,2]<-weight
  x[,3]<-pos
  x[,4]<-groupID
  res<-list(x=x, center=center)
}


optimizeTab<-function(x, tables, gamma)
{
  d <- dim(tables)[1]
  value<-rep(0,d)
  for(i in 1:d)
  {
    t<-tables[i,]
    s <- 0
    for (j in 1:max(t))
    {
      allV<-x[t==j]
      if (length(allV)>0)
      {
        s <- s + gamma
        if (length(allV)>1){
          s<-s + sum((allV-mean(allV))* (allV-mean(allV)))
        }
      }
    }
    value[i]<-s
  }
  ind<-which(value==min(value))[1]
  ind
}


optimizeTabMedian<-function(x, tables)
{
  d <- dim(tables)[1]
  value<-rep(0,d)
  for(i in 1:d)
  {
    t<-tables[i,]
    s <- 0
    cluster<-unique(t)
    value[i]<-length(cluster)
  }
  ind<-which(value==median(value))[1]
  ind  
}


# getDunnIndex<-function(x, result, ind)
# {
#   index<-dunn(clusters=result$tables[ind,], Data = x, method = "euclidean")
#   return (index)
# }


optimizeSil<-function(x, tables)
{
  d <- dim(tables)[1]
  value<-rep(0,d)
  #browser()
  for(i in 1:d)
  {
    t<-tables[i,]
    #value[i]<-dunn(clusters=tables[i,], Data = x, method = "euclidean")
    diss <- as.matrix(dist(x))
    ret <- silhouette(t, dmatrix=diss)

    if(!is.na(ret[1]))
    {
      value[i] <- mean(ret[,3])
    }else{
      value[i] <- -0.1
    }
  }
  ind<-which(value==max(value))[1]
  list(ind=ind, sil=max(value))  
}




