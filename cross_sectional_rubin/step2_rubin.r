## Rubin's Rules for combining results across multiple imputed datasets
## This script loads pre-computed network results and applies Rubin's rules.

# 1. LOAD REQUIRED PACKAGES
library(qgraph)
library(NetworkComparisonTest)
library(openxlsx)
library(officer)
library(flextable)
library(writexl)

# 2. SET PATHS
input_file  <- "results_per_dataset.RData"
output_dir  <- "rubin_output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

datasets_to_process <- c(
  "baseline_all", "week8_all",
  "baseline_rem", "week8_rem",
  "baseline_nonrem", "week8_nonrem"
)

# 3. LOAD PRE-COMPUTED RESULTS
load(input_file)

# 4. HELPER FUNCTIONS FOR DATA CLEANING
extract_clean_values <- function(df_list, ref_names) {
  lapply(df_list, function(df) {
    if (is.null(df) || !("node" %in% names(df)) || !("value" %in% names(df))) {
      return(rep(NA_real_, length(ref_names)))
    }
    df$node <- as.character(df$node)
    ordered_value <- df$value[match(ref_names, df$node)]
    return(ordered_value)
  })
}

extract_predictability_values <- function(df_list, ref_names) {
  lapply(df_list, function(df) {
    if (is.null(df) || !("Variable" %in% names(df)) || !("R2" %in% names(df))) {
      return(rep(NA_real_, length(ref_names)))
    }
    df$Variable <- as.character(df$Variable)
    ordered_r2 <- df$R2[match(ref_names, df$Variable)]
    return(ordered_r2)
  })
}

# 5. RUBIN'S MERGE FUNCTION
merge_with_rubin <- function(results_list, conf_level = 0.95) {
  n_results <- length(results_list)
  first_mat <- results_list[[1]]$edge_matrix
  n_nodes <- nrow(first_mat)
  colnames_abb <- colnames(first_mat)
  z <- qnorm((1 + conf_level) / 2)

  # Initialize storage
  edge_matrices <- array(0, dim = c(n_nodes, n_nodes, n_results))
  ei_list <- strength_list <- predictability_list <- list()
  cs_values <- numeric(n_results)

  for (i in seq_along(results_list)) {
    current <- results_list[[i]]
    edge_matrices[, , i] <- current$edge_matrix
    ei_list[[i]] <- current$centrality
    strength_list[[i]] <- current$strength
    predictability_list[[i]] <- current$predictability
    cs_values[i] <- current$cs
  }

  ref_names <- colnames_abb
  ei_cleaned <- extract_clean_values(ei_list, ref_names)
  strength_cleaned <- extract_clean_values(strength_list, ref_names)
  predict_cleaned <- extract_predictability_values(predictability_list, ref_names)

  # Rubin merge for edge matrices
  rubin_for_edges <- function(edge_mats, z = 1.96) {
    M <- dim(edge_mats)[3]
    mat <- do.call(cbind, lapply(seq_len(M), function(i) as.vector(edge_mats[, , i])))
    theta_bar <- rowMeans(mat)
    B <- apply(mat, 1, var)
    T <- (1 + 1 / M) * B
    se <- sqrt(T)
    ci_low <- theta_bar - z * se
    ci_high <- theta_bar + z * se
    significant <- (ci_low > 0 | ci_high < 0)
    data.frame(
      Estimate = theta_bar,
      SE = se,
      CI_lower = ci_low,
      CI_upper = ci_high,
      Significant = significant
    )
  }

  edge_df <- rubin_for_edges(edge_matrices)
  edge_df$NodeFrom <- rep(colnames_abb, each = n_nodes)
  edge_df$NodeTo <- rep(colnames_abb, times = n_nodes)
  edge_df <- edge_df[, !(names(edge_df) %in% c("From", "To"))]

  # Rubin merge for 1D metrics
  merge_1d <- function(values_list) {
    mat <- do.call(cbind, values_list)
    theta_bar <- rowMeans(mat, na.rm = TRUE)
    B <- apply(mat, 1, var, na.rm = TRUE)
    T <- (1 + 1 / ncol(mat)) * B
    se <- sqrt(T)
    ci_low <- theta_bar - z * se
    ci_high <- theta_bar + z * se
    significant <- (ci_low > 0 | ci_high < 0)
    data.frame(
      Node = ref_names,
      Estimate = theta_bar,
      SE = se,
      CI_lower = ci_low,
      CI_upper = ci_high,
      Significant = significant
    )
  }

  ei_df <- merge_1d(ei_cleaned)
  strength_df <- merge_1d(strength_cleaned)
  predict_df <- merge_1d(predict_cleaned)

  # CS-coefficient merge
  cs_theta <- mean(cs_values, na.rm = TRUE)
  cs_se <- sd(cs_values, na.rm = TRUE) / sqrt(sum(!is.na(cs_values)))
  cs_ci_low <- cs_theta - z * cs_se
  cs_ci_high <- cs_theta + z * cs_se
  cs_significant <- (cs_ci_low > 0 | cs_ci_high < 0)

  cs_df <- data.frame(
    Metric = "CS-Coefficient",
    Estimate = cs_theta,
    SE = cs_se,
    CI_lower = cs_ci_low,
    CI_upper = cs_ci_high,
    Significant = cs_significant
  )

  list(
    EdgeMatrix = edge_df,
    ExpectedInfluence = ei_df,
    Strength = strength_df,
    Predictability = predict_df,
    CStability = cs_df
  )
}

# 6. MAIN LOOP: PROCESS EACH DATASET TYPE
for (dataset_name in datasets_to_process) {
  
  cat("\n Processing dataset:", dataset_name, "\n")
  
  # Extract valid results for this dataset type
  valid_results <- lapply(names(results_per_dataset), function(sheet_name) {
    sheet <- results_per_dataset[[sheet_name]]
    if (!is.null(sheet) && !is.null(sheet[[dataset_name]]) && !is.null(sheet[[dataset_name]]$edge_matrix)) {
      return(sheet[[dataset_name]])
    } else {
      return(NULL)
    }
  })
  
  valid_results <- Filter(Negate(is.null), valid_results)
  
  if (length(valid_results) == 0) {
    cat("Dataset '", dataset_name, "' is missing or invalid. Skipping...\n")
    next
  }
  
  # Align node order
  ref_names <- colnames(valid_results[[1]]$edge_matrix)
  if (is.null(ref_names)) stop("First dataset missing column names")
  
  aligned_results <- lapply(valid_results, function(res) {
    mat <- res$edge_matrix
    if (!identical(colnames(mat), ref_names)) {
      if (all(ref_names %in% colnames(mat))) {
        res$edge_matrix <- mat[ref_names, ref_names]
        return(res)
      } else {
        warning(" Missing nodes:", paste(setdiff(ref_names, colnames(mat)), collapse = ", "))
        return(NULL)
      }
    } else {
      return(res)
    }
  })
  
  aligned_results <- Filter(Negate(is.null), aligned_results)
  
  if (length(aligned_results) < 2) {
    cat("At least two consistent datasets required for Rubin merging '", dataset_name, "'. Skipping...\n")
    next
  }
  
  # Perform Rubin merge
  merged_results <- merge_with_rubin(aligned_results)
  
  # Save results
  output_file <- file.path(output_dir, paste0("rubin_merged_", dataset_name, ".xlsx"))
  write_xlsx(merged_results, path = output_file)
  
  cat("Dataset '", dataset_name, "' merged and saved to:", output_file, "\n")
}