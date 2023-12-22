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

# Similarity Matrix: element (i,j) is the power parameter between baskets i and j
localPP(nDat = c(40,40,40,40,40), yDat = c(3,10,12,18,20), 
        be.a0 = rep(0.5, 5), be.b0 = rep(0.5, 5), a = 1, delta = 1)$sm
