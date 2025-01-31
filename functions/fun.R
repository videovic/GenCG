###--------------------------- Functions require igraph, other packages

run_experiments <- function(experiments) {
  for (experiment_name in names(experiments)) {
    experiment <- experiments[[experiment_name]]
    func <- experiment$func
    params <- experiment$params
    message("Running ", experiment_name)
    do.call(func, params)
  }
}

get_mixing_matrix <- function(g) {
  degrees <- degree(g)
  degrees <- degrees[degrees>0]  # remove isolates
  excess_degrees <- degrees - 1
  
  max_excess_degree <- max(excess_degrees)
  mixing_matrix <- matrix(0, nrow = max_excess_degree + 1, ncol = max_excess_degree + 1)
  
  for (e in E(g)) {
    u <- ends(g, e)[1]
    v <- ends(g, e)[2]
    
    u_excess_degree <- excess_degrees[u] + 1 # R starts at 0
    v_excess_degree <- excess_degrees[v] + 1
    
    mixing_matrix[u_excess_degree, v_excess_degree] <- mixing_matrix[u_excess_degree, v_excess_degree] + 1
    mixing_matrix[v_excess_degree, u_excess_degree] <- mixing_matrix[v_excess_degree, u_excess_degree] + 1
  }
  
  total_edges <- sum(mixing_matrix)
  mixing_matrix <- mixing_matrix / total_edges
  
  return(mixing_matrix)
}


generate_theta_combinations <- function(funs, thetas) {
  param_combinations <- expand.grid(
    lapply(thetas, function(theta) if (is.null(theta)) NULL else theta), # covert to NA -> NULL
    KEEP.OUT.ATTRS = FALSE
  )
  colnames(param_combinations) <- names(funs)
  
  param_combinations <- lapply(1:nrow(param_combinations), function(i) {
    as.list(param_combinations[i, ])
  })
  
  return(param_combinations)
}
