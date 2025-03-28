lagged_analysis <- function(funs, n, thetas, lags,n_reps, ttl,  filename){

    all_assortativity_values <- list()
    
    for (i in seq_along(lags)) {
      lag <- lags[i]
      theta0s <- thetas - lag
      
      result_matrix <- matrix(NA, nrow = length(thetas), ncol = n_reps)
      
      for (j in seq_along(thetas)) {
        result_matrix[, j] <- replicate(n_reps, assortativity.degree(graphon_tensor(funs, n, list(thetas[j], theta0s[j]))))
      }
      
      all_assortativity_values[[i]] <- result_matrix
    }
    
    write.table(all_assortativity_values, "all_assortativity_values.txt")
    write.csv(all_assortativity_values, "all_assortativity_values.csv")
    
    simulation_df <- do.call(rbind, lapply(seq_along(all_assortativity_values), function(i) {
      result_matrix <- all_assortativity_values[[i]]
      
      df_list <- lapply(1:nrow(result_matrix), function(k) {
        data.frame(
          x = rep(thetas[k], n_reps),
          y = result_matrix[, k],
          lag = lags[i]
        )
      })
      
      do.call(rbind, df_list)
    }))
    
    colnames(simulation_df) <- c("x", "y", "lag")
    
    summary_stats <- simulation_df %>%
      group_by(x, lag) %>%
      summarise(
        average = mean(y),
        min = min(y),
        max = max(y),
        .groups = "drop"
      )
    
    # Generate colours based on the number of lags
    num_col <- length(lags)
    cols <- brewer.pal(min(num_col, length(lags)), "Set1")
    
    summary_stats$color <- rep(cols, length.out = nrow(summary_stats))
    
    write.table(summary_stats, "summary_stats.txt")
    simulation_all <- summary_stats %>% filter(lag %in% c(0,2,4)) %>% filter( x > -6)
    
    p <- ggplot(simulation_all) +
      geom_ribbon(aes(x = x, ymin = min, ymax = max, fill = as.factor(lag)), alpha = 0.3) +
      geom_line(aes(x = x, y = average, color = as.factor(lag)), linewidth = 1.2) +
      geom_point(aes(x = x, y = average, color = as.factor(lag)), size = 3) +
      geom_point(aes(x = x, y = average, color = as.factor(lag), shape = as.factor(lag)), size = 3) +
      labs(
        title = ttl,
        x = expression(theta),
        y = "r",
        color = "Lag",
        shape = "Lag",
        fill = "Lag"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        legend.position = "right"
      ) +
      scale_x_continuous(
        breaks = seq(floor(min(thetas)), ceiling(max(thetas)), by=1)
      ) +
      scale_y_continuous(
        breaks = seq(floor(min(simulation_all$min)), ceiling(max(simulation_all$max)), by = 0.02),
        labels = scales::label_number(accuracy = 0.01)
      ) +
      scale_color_manual(values = cols) +
      scale_fill_manual(values = cols)
    
    p
    
    # Save the plot
    ggsave(filename, p, device = "pdf")
    return(p)

}

lagged_analysis(funs=list(c_frank, c_frank), 
                          n=1000, 
                          thetas = c(-10:-1), 
                          lags = c(0:5),
                          n_reps=10, 
                          ttl = TeX("$W_F \\otimes W_F$ Lagged Assortativity"),
                          filename="outputs/frank-lagged-disassort-02.pdf")

## ---------------------------------------------------------------------------
# Zero-assortativity graphs

get_P3 <- function(g){
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
    
  A_cubed <- A %*% A %*% A
  C3e <- sum(diag(A_cubed))/6
  P3e <- sum_value - 3*C3e
  
  return(P3e)
}

get_P2 <-function(g){
  x<- sum(degree(g)*(degree(g)-1))/2  
  return(x)
}

comparison_analysis_r_graphs <- function(funs, n, theta, n_reps, ttl, labels, ylabs,
                                         filename1, 
                                         filename2){
  
  # sim1 <- matrix(0, nrow = n_reps, ncol = length(theta))
  # sim2 <- sim1
  # 
  # for (j in 1:length(theta)) {
  #   sim1[, j] <- sapply(1:n_reps, function(i) {
  #     assortativity.degree(graphon_simulator(funs[[1]], n, theta[j]))
  #   })
  #   sim2[, j] <- sapply(1:n_reps, function(i) {
  #     assortativity.degree(graphon_simulator(funs[[2]], n, theta[j]))
  #   })
  # }
  
  sim_results <- array(list(), dim = c(n_reps, length(theta), 2))
  for (j in 1:length(theta)) {
    sim_results[, j, 1] <- lapply(1:n_reps, function(i) {
      graphon_simulator(funs[[1]], n, theta[j])
    })
    sim_results[, j, 2] <- lapply(1:n_reps, function(i) {
      graphon_simulator(funs[[2]], n, theta[j])
    })
  }

  assort_metrics <- array(NA, dim = c(n_reps, length(theta), 2))
  for (j in 1:length(theta)) {
    for ( i in 1:n_reps){
      assort_metrics[i, j, 1] <- assortativity.degree(sim_results[[i, j, 1]])
      assort_metrics[i, j, 2] <- assortativity.degree(sim_results[[i, j, 2]])
    }
  }
  
  # assort_stats <- apply(assort_metrics, c(2, 3), function(x) c(min = min(x), mean = mean(x), max = max(x)))
  # assort_stats_matrix <- aperm(assort_stats, c(2, 1, 3))
  
  # P3_est <- array(NA, dim = c(n_reps, length(theta), 2))
  # for (j in 1:length(theta)) {
  #   for ( i in 1:n_reps){
  #     P3_est[i, j, 1] <- get_P3(sim_results[[i, j, 1]])
  #     P3_est[i, j, 2] <- get_P3(sim_results[[i, j, 2]])
  #   }
  # }
  
  # P2_est <- array(NA, dim = c(n_reps, length(theta), 2))
  # for (j in 1:length(theta)) {
  #   for ( i in 1:n_reps){
  #     P2_est[i, j, 1] <- get_P2(sim_results[[i, j, 1]])
  #     P2_est[i, j, 2] <- get_P2(sim_results[[i, j, 2]])
  #   }
  # }
  
  P3_P2_est <- array(NA, dim = c(n_reps, length(theta), 2))
  for (j in 1:length(theta)) {
    for ( i in 1:n_reps){
      P3_P2_est[i, j, 1] <- get_P3(sim_results[[i, j, 1]])/get_P2(sim_results[[i, j, 1]])
      P3_P2_est[i, j, 2] <- get_P3(sim_results[[i, j, 1]])/get_P2(sim_results[[i, j, 2]])
    }
  }
  
  P2_P1_est <- array(NA, dim = c(n_reps, length(theta), 2))
  for (j in 1:length(theta)) {
    for ( i in 1:n_reps){
      P2_P1_est[i, j, 1] <- get_P2(sim_results[[i, j, 1]])/(0.5*sum(degree(sim_results[[i, j, 1]])))
      P2_P1_est[i, j, 2] <- get_P2(sim_results[[i, j, 1]])/(0.5*sum(degree(sim_results[[i, j, 2]])))
    }
  }

  graph_plotter<- function(data, labels, ttl, ylabs, filename){
    
    stats_data <- apply(data, c(2, 3), function(x) c(min = min(x), mean = mean(x), max = max(x)))
    stats_data_matrix <- aperm(stats_data, c(2, 1, 3))
    
      stats_sum_1 <- data.frame(
        x = theta,
        min = stats_data_matrix[, "min",1],
        mean = stats_data_matrix[, "mean",1],
        max = stats_data_matrix[, "max",1]
      )
      
      stats_sum_2<- data.frame(
        x = theta,
        min = stats_data_matrix[, "min",2],
        mean = stats_data_matrix[, "mean",2],
        max = stats_data_matrix[, "max",2]
      )
      
      stats_sum_1$function_label <- labels[[1]]
      stats_sum_2$function_label <- labels[[2]]
      
      plot_data <- rbind(stats_sum_1, stats_sum_2)
      
      p <- ggplot(plot_data) +
        
        geom_ribbon(aes(x = x, ymin = min, ymax = max, fill = function_label), alpha = 0.3) +
        geom_line(aes(x = x, y = mean, color = function_label, linetype = function_label), linewidth = 1.2) +
        
        labs(title = ttl1,
             x = expression(theta),
             y = ylab) +
        
        # scale_fill_manual(values = c("Gumbel" = "lightblue", "Joe" = "lightcoral")) +
        # scale_color_manual(values = c("Gumbel" = "blue", "Joe" = "darkred")) +
        # scale_linetype_manual(values = c("Gumbel" = "dashed", "Joe" = "dashed")) +
        
        scale_fill_manual(values = setNames(c("lightblue", "lightcoral"), labels)) +
        scale_color_manual(values = setNames(c("blue", "darkred"), labels)) +
        scale_linetype_manual(values = setNames(c("dashed", "dashed"), labels)) +
        
        theme_bw() +
        theme(
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 12),
          legend.position = c(0.95, 0.05),  
          legend.direction = "horizontal",  
          legend.title = element_blank(),   
          legend.spacing.x = unit(0.3, "cm"),
          legend.justification = c("right", "bottom")  
        ) +
        
        scale_x_continuous(
          breaks = seq(floor(min(plot_data$x)), ceiling(max(plot_data$x)), by = 1)  # Only show integer breaks
        ) +
        
        guides(
          fill = guide_legend(title = NULL),   
          color = guide_legend(title = NULL),  
          linetype = guide_legend(title = NULL) 
        )
      
      ggsave(filename, p, device = "pdf")
      return(p)
  }
  
  graph_plotter(assort_metrics, labels, ttl1, ylabs[[1]], filename1)
  #graph_plotter(P3_P2_est, labels, ttl2, ylabs[[2]], filename2)
  #graph_plotter(P2_P1_est, labels, ttl2, ylabs[[2]], filename2)
}

comparison_analysis_r_graphs(funs=list(c_gumbel, c_joe), 
                             n=1000, 
                             theta = c(1:5), 
                             n_reps=10, 
                             ttl1 = TeX("$W_G$, $W_J$ Assortativity Coefficient"),
                             ttl2 = TeX("$W_G$, $W_J$ P3"),
                             labels = list("Gumbel", "Joe"),
                             ylabs = list("r", TeX("|P_2| / |P_1|")),
                             filename1="outputs/zero-assort-01.pdf",
                             filename2="outputs/zero-assort-02.pdf")
