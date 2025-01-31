## ---------------------------------------------------------------------------
# Helpful Plotting Funs
analysis_fan_r<- function(funs, n, thetas, n_reps, ttl, filename){

  simulations <- matrix(0, nrow = n_reps, ncol = length(thetas))
  
  for (j in 1:length(thetas)) {
    simulations[, j] <- sapply(1:n_reps, function(i) {
      assortativity.degree(graphon_simulator(funs, n, thetas[j]))
    })
  }
  average <- colMeans(simulations)
  min_vals <- apply(simulations, 2, min)
  max_vals <- apply(simulations, 2, max)
  
  flat_thetas <- unlist(thetas, use.names = FALSE)
  flat_thetas <- flat_thetas[!is.na(flat_thetas)]
  
  plot_data <- data.frame(
    x = rep(flat_thetas, each = n_reps),
    y = as.vector(simulations),
    rep = rep(1:n_reps, times = length(flat_thetas))
  )
  
  summary_data <- data.frame(
    x = flat_thetas,
    average = average,
    min = min_vals,
    max = max_vals
  )
  
  p <- ggplot() +
    
    # shaded range (fan graph)
    geom_ribbon(data = summary_data, aes(x = x, ymin = min, ymax = max), fill = "lightblue", alpha = 0.5) +
    
    # average line
    geom_line(data = summary_data, aes(x = x, y = average), color = "blue", linewidth = 1.2) +
    
    labs(title = ttl,
         x = expression(theta),
         y = "r") +
    theme_bw() +  
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      legend.position = "none"  
    )+
    scale_x_continuous(
      breaks = seq(floor(min(summary_data$x)), ceiling(max(summary_data$x)), by = 1)  # Only show integer breaks
    )
  
  ggsave(filename, p, device = "pdf")
  
  return(p)
  
}

analysis_fan_tensor_r<- function(funs, n, theta1, theta2s, n_reps, ttl, filename){
  
  simulations <- matrix(0, nrow = n_reps, ncol = length(theta2s))
  
  for (j in 1:length(theta2s)) {
    simulations[, j] <- sapply(1:n_reps, function(i) {
      assortativity.degree(graphon_tensor(funs, n,list(theta1, theta2s[j])))
    })
  }
  
  average <- colMeans(simulations)
  min_vals <- apply(simulations, 2, min)
  max_vals <- apply(simulations, 2, max)
  
  flat_thetas <- unlist(theta2s, use.names = FALSE)
  flat_thetas <- flat_thetas[!is.na(flat_thetas)]
  
  plot_data <- data.frame(
    x = rep(flat_thetas, each = n_reps),
    y = as.vector(simulations),
    rep = rep(1:n_reps, times = length(flat_thetas))
  )
  
  summary_data <- data.frame(
    x = flat_thetas,
    average = average,
    min = min_vals,
    max = max_vals
  )
  
  p <- ggplot() +
    
    # shaded range (fan graph)
    geom_ribbon(data = summary_data, aes(x = x, ymin = min, ymax = max), fill = "lightblue", alpha = 0.5) +
    
    # average line
    geom_line(data = summary_data, aes(x = x, y = average), color = "blue", linewidth = 1.2) +
    
    labs(title = ttl,
         x = expression(theta),
         y = "r") +
    theme_bw() +  
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      legend.position = "none"  
    )+
    scale_x_continuous(
      breaks = seq(floor(min(summary_data$x)), ceiling(max(summary_data$x)), by = 1)  # Only show integer breaks
    )
  
  ggsave(filename, p, device = "pdf")
  
  return(p)
  
}

analysis_fan_rep_tensor_r<- function(funs, n, thetas, n_reps, ttl, filename){
  
  simulations <- matrix(0, nrow = n_reps, ncol = length(thetas))
  
  for (j in 1:length(thetas)) {
    simulations[, j] <- sapply(1:n_reps, function(i) {
      assortativity.degree(graphon_tensor(funs, n,list(thetas[j], thetas[j])))
    })
  }
  
  average <- colMeans(simulations)
  min_vals <- apply(simulations, 2, min)
  max_vals <- apply(simulations, 2, max)
  
  flat_thetas <- unlist(thetas, use.names = FALSE)
  flat_thetas <- flat_thetas[!is.na(flat_thetas)]
  
  plot_data <- data.frame(
    x = rep(flat_thetas, each = n_reps),
    y = as.vector(simulations),
    rep = rep(1:n_reps, times = length(flat_thetas))
  )
  
  summary_data <- data.frame(
    x = flat_thetas,
    average = average,
    min = min_vals,
    max = max_vals
  )
  
  p <- ggplot() +
    
    # shaded range (fan graph)
    geom_ribbon(data = summary_data, aes(x = x, ymin = min, ymax = max), fill = "lightblue", alpha = 0.5) +
    
    # average line
    geom_line(data = summary_data, aes(x = x, y = average), color = "blue", linewidth = 1.2) +
    
    labs(title = ttl,
         x = expression(theta),
         y = "r") +
    theme_bw() +  
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      legend.position = "none"  
    )+
    scale_x_continuous(
      breaks = seq(floor(min(summary_data$x)), ceiling(max(summary_data$x)), by = 1)  # Only show integer breaks
    )
  
  ggsave(filename, p, device = "pdf")
  
  return(p)
  
}
## ---------------------------------------------------------------------------
# 1- No Assortativity

# graph $r$ in a small neighborhood around $r=0$
# Note J,G as \theta ->1, r -> 0


## ---------------------------------------------------------------------------
# 2- Assortativity

#weak
analysis_fan_r(funs=c_gumbel, 
               n=1000, 
               thetas = c(1:10), 
               n_reps=10, 
               ttl= expression( paste(W[F], " Gumbel Copula Assortativity")),
               filename="outputs/gumbel-weakassort-01.pdf")

#wider assortativity range
analysis_fan_tensor_r(funs=list(d_joe, d_gumbel), 
               n=1000, 
               theta1= 2,
               theta2s = c(1:10), 
               n_reps=10, 
               ttl = TeX("$\\tilde{W}_F \\otimes \\tilde{W}_G$ Assortativity"),
               filename="outputs/gumbeljoe-highassort-01.pdf")

## ---------------------------------------------------------------------------
# 3- Disassortativity

analysis_fan_rep_tensor_r(funs=list(c_frank, c_frank), 
                      n=1000, 
                      thetas = c(-10:-1), 
                      n_reps=10, 
                      ttl = TeX("$W_F \\otimes W_F$ Assortativity"),
                      filename="outputs/frank-disassort-01.pdf")

