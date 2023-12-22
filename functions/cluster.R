library(partitions)
# These functions are required for fitting local MEM (Liu et al. 2022)
# R: the number of trials to consider
# max_cl: maximum number of blocks
# same results as previous version
get.part <- function(R, max_cl){
  ## generate all the possible partitions 
  ## and store them in a matrix 
  part_mat <- t(setparts(R))
  part_mat <- part_mat[apply(part_mat, 1, function(x){
    length(unique(x))<=max_cl
  }), ]
  part_mat <- data.frame(part_mat)
  names(part_mat) <- LETTERS[1:R]
  #part_mat[-1, ]
  part_mat
}

update.part <- function(x, n, prior_part, part, a0 = 0.5, b0 = 0.5){
  R <- length(x)
  K <- nrow(part)
  
  p <- foreach(k = 1:K,.combine = "c")%do%{
    grp <- unlist(part[k, ])
    S <- aggregate(x, by=list(grp), sum)$x
    N <- aggregate(n, by=list(grp), sum)$x
    #calculate marginal probs m(s_j)
    prod((beta(a0+S, b0+N-S)/beta(a0,b0)))*prior_part[k]
  }
  #marginal probs m(s_j)
  mp <- sum(p)
  ## calculate posterior prob of each grouping structure
  p/mp
}

