#fun.R

## Generates graphons by a single f
## Lovasz, 2012

graphon_simulator <-function(f,n,theta0){ 
  p <- A <- matrix(0, n, n)
  u <-runif(n)
  for (i in 1:(n - 1)) {
    for (j in (i+1):n) {
      p[i, j] <- p[j, i] <- f(u[i], u[j], theta0)
      A [i, j] <- A [j, i] <- runif(1) < p[i, j] 
    }
  }
  g <- graph_from_adjacency_matrix(A, "undirected", diag=FALSE)
  return(g)
}

graphon_simulator_fundamental <-function(f,n){ 
  p <- A <- matrix(0, n, n)
  u <-runif(n)
  for (i in 1:(n - 1)) {
    for (j in (i+1):n) {
      p[i, j] <- p[j, i] <- f(u[i], u[j])
      A [i, j] <- A [j, i] <- runif(1) < p[i, j] 
    }
  }
  g <- graph_from_adjacency_matrix(A, "undirected", diag=FALSE)
  return(g)
}

## A function to generate tensor based combinations of graphons
## Lovasz, 2012
graphon_tensor <- function(fs, n, theta0){
  if (!is.list(fs)) fs <- list(fs)
  if (!is.list(theta0)) theta0 <- list(theta0)

  if ( length(fs) != length(theta0)) {
    stop("The lengths of fs and theta0 must be the same.")
  }
  
  # Partitions of uniforms
  num_functions <- length(fs)
  u_list <- lapply(seq_along(fs), function(k) runif(n))

  p <- matrix(1, nrow = length(u_list[[1]]), ncol = length(u_list[[1]]))
  # set to 1, so that it will multiply something
  
  p_i <- matrix(0, nrow = length(u_list[[1]]), ncol = length(u_list[[1]]))
  
  for (k in seq_along(fs)) {
    f <- fs[[k]]
    theta00 <- theta0[[k]]
    u_k <- u_list[[k]]
    
    if(is.na(theta00)){
      p_i <- outer(u_k, u_k, Vectorize(function(x, y){ f(x,y)}))
    } else {
      p_i <- outer(u_k, u_k, Vectorize(function(x, y,theta){ f(x,y, theta00)}))
    }
    
    p <- p*p_i# element-wise multiplication
  }

  A <- matrix(0, n, n)
  A[upper.tri(A)] <- runif(sum(upper.tri(A))) < p[upper.tri(p)]  
  A <- A + t(A) 
  
  g <- graph_from_adjacency_matrix(A, mode = "undirected", diag=FALSE)
  return(g)
}





