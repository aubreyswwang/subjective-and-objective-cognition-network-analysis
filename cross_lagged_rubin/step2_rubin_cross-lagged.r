## Rubin's Rules for combining cross-lagged panel network results across multiple imputed datasets

# 1. SET CRAN MIRROR
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

# 2. LOAD REQUIRED PACKAGES
pkgs <- c("mice", "openxlsx", "dplyr", "tidyr", "broom", "writexl")
invisible(lapply(pkgs, function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}))

# 3. SET PATHS (use relative paths)
input_file  <- "cross_lagged_network_results.RData"
output_file <- "rubin_results_cross_lagged.xlsx"

# 4. LOAD DATA
if (!file.exists(input_file)) {
  stop(" Input file not found: ", input_file)
}
load(input_file)

# Validate data structure
required_lists <- c("edges_list", "outEI_list", "inEI_list", "cs_list")
missing_lists <- required_lists[!sapply(required_lists, exists)]
if (length(missing_lists) > 0) {
  stop(" Missing required data lists: ", paste(missing_lists, collapse = ", "))
}

cat(" Data structure validated.\n")
cat(" Available imputed datasets:", paste(names(edges_list), collapse = ", "), "\n")

# 5. RUBIN MERGE FUNCTION
rubin_merge <- function(values) {
  values <- as.numeric(na.omit(values))
  if (length(values) < 2) {
    return(list(Estimate = NA, SE = NA, CI_lower = NA, CI_upper = NA, P_value = NA, Significant = NA))
  }
  
  m <- length(values)
  Q_bar <- mean(values)
  U <- var(values)
  B <- var(values)
  T <- U + (1 + 1/m) * B
  SE <- sqrt(T)
  
  # P-value calculation
  if (SE > 0) {
    t_stat <- Q_bar / SE
    lambda <- (1 + 1/m) * B / T
    df <- (m - 1) * (1 + 1/lambda)^2
    df <- max(df, 1)
    p_value <- 2 * pt(abs(t_stat), df = df, lower.tail = FALSE)
  } else {
    p_value <- NA
  }
  
  # Confidence interval
  t_crit <- qt(0.975, df = df)
  CI_lower <- Q_bar - t_crit * SE
  CI_upper <- Q_bar + t_crit * SE
  Significant <- !between(0, CI_lower, CI_upper)
  
  list(
    Estimate = Q_bar,
    SE = SE,
    CI_lower = CI_lower,
    CI_upper = CI_upper,
    P_value = p_value,
    Significant = Significant
  )
}

# 6. EXTRACT AND ALIGN DATA
safe_extract <- function(df, node_name, column_name) {
  if (is.data.frame(df) && node_name %in% row.names(df)) {
    return(df[node_name, column_name])
  } else {
    return(NA_real_)
  }
}

all_nodes <- unique(unlist(lapply(inEI_list, row.names)))

inEI_values <- lapply(names(inEI_list), function(sheet) {
  sapply(all_nodes, function(node) safe_extract(inEI_list[[sheet]], node, "InEI.value"))
})
outEI_values <- lapply(names(outEI_list), function(sheet) {
  sapply(all_nodes, function(node) safe_extract(outEI_list[[sheet]], node, "OutEI.value"))
})

inEI_values <- t(do.call(rbind, inEI_values))
outEI_values <- t(do.call(rbind, outEI_values))
rownames(inEI_values) <- rownames(outEI_values) <- all_nodes
colnames(inEI_values) <- colnames(outEI_values) <- names(inEI_list)

# Extract edge values and CS coefficients
edge_values <- lapply(names(edges_list), function(sheet) {
  edges_list[[sheet]]$edge.value
})
cs_inEI_values <- lapply(cs_list, function(x) x$Stability_Coefficient[1])
cs_outEI_values <- lapply(cs_list, function(x) x$Stability_Coefficient[2])

# 7. APPLY RUBIN'S RULES
# Edges
edge_merged <- lapply(seq_along(edge_values[[1]]), function(i) {
  values <- sapply(edge_values, `[[`, i)
  rubin_merge(values)
})
edges_summary <- edges_list[[1]][, c("edge.nodeOut", "edge.nodeIn")]
edges_summary$Estimate <- sapply(edge_merged, `[[`, "Estimate")
edges_summary$SE <- sapply(edge_merged, `[[`, "SE")
edges_summary$CI_lower <- sapply(edge_merged, `[[`, "CI_lower")
edges_summary$CI_upper <- sapply(edge_merged, `[[`, "CI_upper")
edges_summary$P_value <- sapply(edge_merged, `[[`, "P_value")
edges_summary$Significant <- sapply(edge_merged, `[[`, "Significant")

# InEI
inEI_merged <- apply(inEI_values, 1, rubin_merge)
inEI_summary <- data.frame(Node = rownames(inEI_values))
inEI_summary$Estimate <- sapply(inEI_merged, `[[`, "Estimate")
inEI_summary$SE <- sapply(inEI_merged, `[[`, "SE")
inEI_summary$CI_lower <- sapply(inEI_merged, `[[`, "CI_lower")
inEI_summary$CI_upper <- sapply(inEI_merged, `[[`, "CI_upper")
inEI_summary$P_value <- sapply(inEI_merged, `[[`, "P_value")
inEI_summary$Significant <- sapply(inEI_merged, `[[`, "Significant")

# OutEI
outEI_merged <- apply(outEI_values, 1, rubin_merge)
outEI_summary <- data.frame(Node = rownames(outEI_values))
outEI_summary$Estimate <- sapply(outEI_merged, `[[`, "Estimate")
outEI_summary$SE <- sapply(outEI_merged, `[[`, "SE")
outEI_summary$CI_lower <- sapply(outEI_merged, `[[`, "CI_lower")
outEI_summary$CI_upper <- sapply(outEI_merged, `[[`, "CI_upper")
outEI_summary$P_value <- sapply(outEI_merged, `[[`, "P_value")
outEI_summary$Significant <- sapply(outEI_merged, `[[`, "Significant")

# CS coefficients
cs_inEI_merged <- rubin_merge(unlist(cs_inEI_values))
cs_outEI_merged <- rubin_merge(unlist(cs_outEI_values))
cs_summary <- data.frame(
  Centrality = c("inExpectedInfluence", "outExpectedInfluence"),
  Estimate = c(cs_inEI_merged$Estimate, cs_outEI_merged$Estimate),
  SE = c(cs_inEI_merged$SE, cs_outEI_merged$SE),
  CI_lower = c(cs_inEI_merged$CI_lower, cs_outEI_merged$CI_lower),
  CI_upper = c(cs_inEI_merged$CI_upper, cs_outEI_merged$CI_upper),
  P_value = c(cs_inEI_merged$P_value, cs_outEI_merged$P_value),
  Significant = c(cs_inEI_merged$Significant, cs_outEI_merged$Significant)
)

# 8. SAVE RESULTS
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
write_xlsx(list(
  edges = edges_summary,
  outEI = outEI_summary,
  inEI = inEI_summary,
  cs = cs_summary
), path = output_file)

cat(" Rubin merging completed. Results saved to:", output_file, "\n")