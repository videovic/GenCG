# counts of homomorphism numbers/ densities and useful isomorphism numbers / densities
# requires library(cubature) 

####################################################################
# gives theoretical densities (estimated) via integral functions
# for some homomorphism densities
compute_S_k <- function(k, fun, theta = NULL) {
  # star on k+1 nodes and k edges
  
  inner_function <- function(x, y) {
    if (is.null(theta)) {
      fun(x, y)^k
    } else {
      fun(x, y, theta)^k
    }
  }
  outer_integral <- function(y) {
    adaptIntegrate(function(x) inner_function(x, y),
                   lowerLimit = 0,
                   upperLimit = 1, tol = 1e-8)$integral
  }
  result <- adaptIntegrate(function(y) outer_integral(y),
                           lowerLimit = 0,
                           upperLimit = 1, tol = 1e-8)$integral
  return(result)
}

compute_Pk <- function(k, fun, theta = NULL) {
  #on k edges and k+1 nodes
  integrand <- function(vars) {
    product <- 1
    for (i in 1:k) {
      if (is.null(theta)) {
        product <- product * fun(vars[i], vars[i + 1])
      } else {
        product <- product * fun(vars[i], vars[i + 1], theta)
      }
    }
    return(product)
  }
  
  result <- adaptIntegrate(integrand,
                           lowerLimit = rep(0, k + 1),
                           upperLimit = rep(1, k + 1),
                           tol = 1e-8)$integral
  
  return(result)
}

compute_C3 <- function(fun, theta = NULL) {
  integrand <- function(vars) {
    x <- vars[1]
    y <- vars[2]
    z <- vars[3]
    if (is.null(theta)) {
      fun(x, y) * fun(y, z) * fun(z, x)
    } else {
      fun(x, y, theta) * fun(y, z, theta) * fun(z, x, theta)
    }  
  }
  result <- adaptIntegrate(integrand,
                           lowerLimit = c(0, 0, 0),
                           upperLimit = c(1, 1, 1), tol = 1e-8)$integral
  return(result)
}

####################################################################
# gives empirical estimates of functionals, based off an observed g

estimate_P1 <- function(g) 0.5*sum(degree(g)) 

# gives injective mappings for a path of length k
# note P2 = S2
estimate_Pk_inj <- function(g,k) sum(choose(degree(g), k))

P3_sum <- function(g){
  A <- as.matrix(as_adjacency_matrix(g))
  degrees <- rowSums(A)
  sum_value <- 0
  
  for (i in 1:nrow(A)) {
    for (j in i:ncol(A)) {
      if (A[i, j] == 1) {
        k_i <- degrees[i]
        k_j <- degrees[j]
        sum_value <- sum_value + (k_i - 1) * (k_j - 1)
      }
    }
  }
  return(sum_value)
}

P3 <- function(g){
  A <- as.matrix(as_adjacency_matrix(g))
  A_cubed <- A %*% A %*% A
  C3 <- sum(diag(A_cubed))/6
  return(P3_sum(g) - 3*C3)
}


estimate_S_13_hom  <- function(g) sum((degree(g)^3))

estimate_tSk <- function(g,k)  sum((degree(g)^k))/(gorder(g)^(k+1))

clustering_coeff <- function(g) transitivity(g, type = "global")


