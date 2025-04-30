# Load required packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  dplyr,           # For data manipulation
  tidyr,           # For data tidying
  ggplot2,         # For plotting
  patchwork,       # For combining plots
  pryr,            # For memory_size
  future,          # For parallel processing
  future.apply,    # For parallel apply functions
  scales,          # For nicer plot scales
  gaudi
)

# Configure parallel processing
future::plan(future::multisession)

# Define benchmark parameters
benchmark_config <- expand.grid(
  n_samples = c(50, 100, 200, 500, 1000),       # Number of samples
  n_features = c(1000, 5000, 10000, 25000),     # Number of features per omics
  n_omics = c(2, 3, 5),                         # Number of omics datasets
  repetitions = 3                               # Number of repetitions per config
)

# Function to generate synthetic multi-omics datasets
generate_synthetic_omics <- function(n_samples, n_features_per_omic, n_omics, 
                                     correlation = 0.5, signal_strength = 0.8) {
  
  # Create a list to store the omics datasets
  omics_list <- list()
  
  # Generate a shared latent structure
  latent_dim <- 5
  latent_factors <- matrix(rnorm(n_samples * latent_dim), nrow = n_samples)
  
  # Generate omics-specific datasets with some correlation to the latent factors
  for (i in 1:n_omics) {
    # Create feature profiles with some biological relevance
    # Each feature has some correlation with the latent factors
    feature_loadings <- matrix(rnorm(n_features_per_omic * latent_dim), 
                               nrow = n_features_per_omic) * signal_strength
    
    # Generate the base data from the latent structure
    base_signal <- latent_factors %*% t(feature_loadings)
    
    # Add noise to make it realistic
    noise <- matrix(rnorm(n_samples * n_features_per_omic), 
                    nrow = n_samples) * (1 - correlation)
    
    # Combine signal and noise
    omic_data <- base_signal + noise
    
    # Add feature and sample names
    colnames(omic_data) <- paste0("feature", i, "_", 1:n_features_per_omic)
    rownames(omic_data) <- paste0("sample_", 1:n_samples)
    
    omics_list[[i]] <- as.data.frame(omic_data)
  }
  
  return(omics_list)
}

# Function to benchmark GAUDI on a given configuration
benchmark_gaudi <- function(n_samples, n_features, n_omics, seed = 42) {
  set.seed(seed)
  
  # Generate synthetic data
  omics_list <- generate_synthetic_omics(
    n_samples = n_samples,
    n_features_per_omic = n_features,
    n_omics = n_omics
  )
  
  # Configure GAUDI parameters
  umap_params <- list(n_neighbors = 10, n_components = 3, min_dist = 0.1)
  umap_params_conc <- list(n_neighbors = 10, n_components = 2, min_dist = 0.1)
  min_pts <- max(floor(0.05 * n_samples), 5)
  
  # Measure memory before
  mem_before <- pryr::mem_used()
  
  # Measure execution time
  start_time <- Sys.time()
  
  tryCatch({
    # Run GAUDI
    result <- gaudi(
      omics = omics_list,
      umap_params = umap_params,
      umap_params_conc = umap_params_conc,
      min_pts = min_pts,
      compute_features = TRUE,
      reassign_cluster_zero = TRUE,
      samples_in_rows = TRUE
    )
    
    # Measure execution time
    end_time <- Sys.time()
    
    # Measure memory after
    mem_after <- pryr::mem_used()
    
    # Calculate memory usage
    mem_used <- mem_after - mem_before
    
    # Calculate execution time
    exec_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    # Return results
    return(list(
      success = TRUE,
      samples = n_samples,
      features = n_features,
      omics = n_omics,
      execution_time = exec_time,
      memory_used = mem_used,
      n_clusters = length(unique(result@clusters)),
      silhouette = result@silhouette_score
    ))
  }, error = function(e) {
    end_time <- Sys.time()
    exec_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    return(list(
      success = FALSE,
      samples = n_samples,
      features = n_features,
      omics = n_omics,
      execution_time = exec_time,
      memory_used = NA,
      n_clusters = NA,
      silhouette = NA,
      error = as.character(e)
    ))
  })
}

# Run all benchmarks
cat("Starting benchmark with", nrow(benchmark_config), "configurations...\n")

# Create a results dataframe
results <- data.frame()

# Loop through configurations
for (i in 1:nrow(benchmark_config)) {
  config <- benchmark_config[i, ]
  cat(sprintf("Running benchmark %d/%d: %d samples, %d features, %d omics, repetition %d\n", 
              i, nrow(benchmark_config), 
              config$n_samples, config$n_features, config$n_omics, config$repetitions))
  
  # Run benchmark with current configuration
  benchmark_result <- benchmark_gaudi(
    n_samples = config$n_samples,
    n_features = config$n_features, 
    n_omics = config$n_omics,
    seed = i  # Different seed for each run
  )
  
  # Add configuration to results
  result_row <- data.frame(
    config_id = i,
    n_samples = config$n_samples,
    n_features = config$n_features,
    n_omics = config$n_omics,
    repetition = config$repetitions,
    success = benchmark_result$success,
    execution_time = benchmark_result$execution_time,
    memory_used = as.numeric(benchmark_result$memory_used),
    n_clusters = benchmark_result$n_clusters,
    silhouette = benchmark_result$silhouette
  )
  
  results <- rbind(results, result_row)
  
  # Write intermediate results to file
  write.csv(results, "./results_scalability_benchmark/gaudi_benchmark_results.csv", row.names = FALSE)
  
  # Sleep between runs
  Sys.sleep(1)
}

# Calculate summary statistics
summary_results <- results %>%
  dplyr::filter(success == TRUE) %>%
  dplyr::group_by(n_samples, n_features, n_omics) %>%
  dplyr::summarize(
    mean_time = mean(execution_time),
    sd_time = sd(execution_time),
    mean_memory = mean(memory_used, na.rm = TRUE),
    sd_memory = sd(memory_used, na.rm = TRUE),
    mean_clusters = mean(n_clusters, na.rm = TRUE),
    mean_silhouette = mean(silhouette, na.rm = TRUE),
    n_runs = n(),
    .groups = "drop"
  )

# Create a summary of failed runs
failed_runs <- results %>%
  filter(success == FALSE) %>%
  select(n_samples, n_features, n_omics, repetition)

# Save results
write.csv(results, "./results_scalability_benchmark/gaudi_benchmark_results.csv", row.names = FALSE)
write.csv(summary_results, "./results_scalability_benchmark/gaudi_benchmark_summary.csv", row.names = FALSE)
if (nrow(failed_runs) > 0) {
  write.csv(failed_runs, "./results_scalability_benchmark/gaudi_benchmark_failures.csv", row.names = FALSE)
}

# ===== Visualization =====

# Time complexity by number of samples
plot_time_by_samples <- ggplot(summary_results, aes(x = n_samples, y = mean_time, color = factor(n_omics))) +
  geom_point(size = 3) +
  geom_line(aes(group = interaction(n_features, n_omics))) +
  # geom_errorbar(aes(ymin = mean_time - sd_time, ymax = mean_time + sd_time), width = 0.1) +
  facet_wrap(~ n_features, scales = "free_y", labeller = labeller(n_features = function(x) paste0("# Features: ", x))) +
  scale_color_viridis_d(name = "Number of Datasets") +
  scale_y_log10(labels = scales::comma) +
  scale_x_log10(labels = scales::comma) +
  labs(
    title = "GAUDI Time Complexity by Number of Samples",
    x = "Number of Samples (log10 scale)",
    y = "Execution Time in Seconds (log10 scale)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom"
  )

# Time complexity by number of features
plot_time_by_features <- ggplot(summary_results, aes(x = n_features, y = mean_time, color = factor(n_omics))) +
  geom_point(size = 3) +
  geom_line(aes(group = interaction(n_samples, n_omics))) +
  # geom_errorbar(aes(ymin = mean_time - sd_time, ymax = mean_time + sd_time), width = 0.1) +
  facet_wrap(~ n_samples, scales = "free_y", labeller = labeller(n_samples = function(x) paste0("# Samples: ", x))) +
  scale_color_viridis_d(name = "Number of Datasets") +
  scale_y_log10(labels = scales::comma) +
  scale_x_log10(labels = scales::comma) +
  labs(
    title = "GAUDI Time Complexity by Number of Features",
    x = "Number of Features (log10 scale)",
    y = "Execution Time in Seconds (log10 scale)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom"
  )

# Memory usage by dataset size
plot_memory <- ggplot(summary_results, aes(x = n_samples * n_features * n_omics, y = mean_memory, color = factor(n_omics))) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed", color = "grey30") +
  scale_color_viridis_d(name = "Number of Datasets") +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10(labels = function(x) scales::comma(x / 1024^2)) +
  labs(
    title = "GAUDI Memory Usage by Dataset Size",
    x = "Total Dataset Size (samples × features × omics, log10 scale)",
    y = "Memory Usage (MB, log10 scale)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom"
  )

# Clustering quality plot
plot_clusters <- ggplot(summary_results, aes(x = n_samples, y = mean_silhouette, color = factor(n_omics))) +
  geom_point(size = 3) +
  geom_line(aes(group = interaction(n_features, n_omics))) +
  facet_wrap(~ n_features, scales = "fixed", labeller = labeller(n_features = function(x) paste0("Features: ", x))) +
  scale_color_viridis_d(name = "Number of Datasets") +
  ylim(0, 1) +
  labs(
    title = "GAUDI Clustering Quality (Silhouette Score)",
    x = "Number of Samples",
    y = "Average Silhouette Score"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom"
  )

# Combine plots
time_plots <- plot_time_by_samples / plot_time_by_features
quality_plots <- plot_memory / plot_clusters

# Save plots
ggsave("./results_scalability_benchmark/gaudi_time_complexity.pdf", time_plots, width = 12, height = 10)
ggsave("./results_scalability_benchmark/gaudi_quality_plots.pdf", quality_plots, width = 12, height = 10)

# Print summary
cat("\n===== Benchmark Summary =====\n")
cat("Configurations tested:", nrow(benchmark_config), "\n")
cat("Successful runs:", sum(results$success), "\n")
cat("Failed runs:", sum(!results$success), "\n")
cat("\nResults saved to CSV files and plots saved as PDFs.\n")

