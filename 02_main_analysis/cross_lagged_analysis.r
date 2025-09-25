## Cross-lagged panel network analysis

# 1. SET WORKING DIRECTORY AND OUTPUT PATH (use relative paths for open-source)
# Note: Assumes this script is run from project root
save_path_cross_lagged <- "results/cross-lagged"
dir.create(save_path_cross_lagged, showWarnings = FALSE, recursive = TRUE)

# 2. LOAD REQUIRED PACKAGES
pkgs <- c("readxl", "dplyr", "haven", "psych", "car", "NetworkComparisonTest",
          "ggplot2", "mgm", "mice", "bootnet", "qgraph", "psychTools", 
          "glmnet", "lavaan", "readr", "EstimateGroupNetwork", "networktools", 
          "NetworkToolbox", "Hmisc", "wCorr", "writexl", "patchwork", "openxlsx",
          "huge", "officer", "flextable", "gridExtra", "tidyr", "forcats", "Rgraphviz")

invisible(lapply(pkgs, function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}))

# 3. LOAD IMPUTED DATA
# Note: Place data file in 'data/' directory
CoV_l_1 <- read_excel("data/mean_imputed_with_single_ipw.xlsx", sheet = "Mean_Imputed")
CoV_l_2 <- CoV_l_1[CoV_l_1$Retained == "retained" & !is.na(CoV_l_1$Retained), ]

# Extract inverse probability weights
weights_global <- CoV_l_2$IPW_w[CoV_l_2$Retained == "retained" & !is.na(CoV_l_2$Retained)]

# 4. PREPARE DATA FOR CROSS-LAGGED NETWORK
w0b <- as.data.frame(CoV_l_2[, 13:22])  # Week 0
w8b <- as.data.frame(CoV_l_2[, 36:45])  # Week 8

# Define node groups and colors
groups_vec <- c(
  rep("Subjective Cognition", 3),
  rep("Objective Cognition", 4),
  rep("Depressive Symptom", 3)
)

group_colors <- c(
  "Subjective Cognition" = "#C8C2E4",
  "Objective Cognition"  = "#9BBBE1",
  "Depressive Symptom"   = "#F09BA0"
)

node_colors <- group_colors[groups_vec]

items_clean <- c(
  "Subjective Attention", "Subjective Memory", "Subjective Executive Functions",
  "Objective Memory", "Objective Executive Functions", "Objective Attention", "Objective Processing Speed",
  "Sad Mood", "General Interest", "Energy"
)

colnames_abb <- c(
  "S-AT", "S-ME", "S-EF",
  "O-ME", "O-EF", "O-AT", "O-PS",
  "SM", "GI", "EN"
)

groups <- split(items_clean, groups_vec)
item_shortnames <- setNames(colnames_abb, items_clean)

# Apply consistent column names
colnames(w0b) <- colnames_abb
colnames(w8b) <- colnames_abb

# Convert to numeric
w0b[, 1:10] <- lapply(w0b[, 1:10], as.numeric)
w8b[, 1:10] <- lapply(w8b[, 1:10], as.numeric)

# 5. CREATE FULL DATASET WITH COVARIATES
w0_8 <- data.frame(w0b, w8b)
w0_8$age <- CoV_l_2$Age
w0_8 <- huge.npn(w0_8)  # Nonparanormal transformation

# Add gender (convert to numeric: assume 1=M, 2=F or similar)
Sex <- as.numeric(CoV_l_2$Sex)
w0_8$gender <- Sex

# 6. ESTIMATE CROSS-LAGGED NETWORK USING GLMNET
numCovar <- 2  # age + gender
k <- 10        # number of nodes per time point

adjMat3 <- matrix(0, nrow = k + numCovar, ncol = k + numCovar)

# Prepare predictor matrix: [T1 nodes + covariates]
x_matrix <- w0_8[, c(1:k, (2*k + 1):(2*k + numCovar))]

# Regress each T2 node on T1 nodes + covariates
for (i in 1:k) {
  set.seed(123)
  fit <- cv.glmnet(
    x           = as.matrix(x_matrix),
    y           = w0_8[, k + i],
    family      = "gaussian",
    alpha       = 1,
    standardize = TRUE,
    weights     = weights_global
  )
  lam <- fit$lambda.min
  # Extract coefficients (skip intercept)
  adjMat3[, i] <- coef(fit, s = lam)[-1]
}

# Create full adjacency matrix with labels
labels <- c(colnames_abb, "gender", "age")
adjMat3_full <- getWmat(adjMat3, nNodes = k + numCovar, labels = labels, directed = TRUE)

# Extract only T1 -> T2 edges (10x10)
adjMat3_T1T2 <- adjMat3_full[1:k, 1:k]

# 7. PLOT CROSS-LAGGED NETWORK (Version 1: with self-loops)
png(file.path(save_path_cross_lagged, "net_with_self_loops.png"),
    width = 12, height = 8, units = "in", res = 300)
par(fig = c(0, 0.65, 0, 1), mar = c(2, 4, 2, 2))
nwc1 <- qgraph(
  adjMat3_T1T2,
  groups = groups_vec,
  legend = FALSE,
  threshold = 0.05,
  labels = colnames_abb,
  color = node_colors,
  directed = TRUE
)

# Add custom legend
par(fig = c(0.65, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0))
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 3), ylim = c(0, 6))
x_left <- 0; x_right <- 2.5; y_pos <- 5.8; item_height <- 0.2
for (group_name in names(groups)) {
  text(x_left + 0.05, y_pos - 0.18, group_name, font = 2, cex = 1.5, adj = c(0, 0.5))
  y_pos <- y_pos - item_height * 2
  for (item in groups[[group_name]]) {
    rect(x_left, y_pos - item_height, x_right, y_pos,
         col = group_colors[group_name], border = "black")
    text(x_left + 0.05, y_pos - item_height / 2,
         paste0(item, " (", item_shortnames[item], ")"),
         adj = c(0, 0.5), cex = 0.9)
    y_pos <- y_pos - item_height
  }
}
dev.off()

# 8. PLOT CROSS-LAGGED NETWORK (Version 2: without self-loops)
adjMat4 <- adjMat3_T1T2
diag(adjMat4) <- 0

png(file.path(save_path_cross_lagged, "net_without_self_loops.png"),
    width = 14, height = 8, units = "in", res = 300)
par(fig = c(0, 0.65, 0, 1), mar = c(2, 4, 2, 2))
nwc2 <- qgraph(
  adjMat4,
  groups = groups_vec,
  legend = FALSE,
  threshold = 0.05,
  labels = colnames_abb,
  color = node_colors,
  directed = TRUE
)

# Same legend
par(fig = c(0.65, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0))
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 3), ylim = c(0, 6))
x_left <- 0; x_right <- 2.5; y_pos <- 5.8; item_height <- 0.2
for (group_name in names(groups)) {
  text(x_left + 0.05, y_pos - 0.18, group_name, font = 2, cex = 1.5, adj = c(0, 0.5))
  y_pos <- y_pos - item_height * 2
  for (item in groups[[group_name]]) {
    rect(x_left, y_pos - item_height, x_right, y_pos,
         col = group_colors[group_name], border = "black")
    text(x_left + 0.05, y_pos - item_height / 2,
         paste0(item, " (", item_shortnames[item], ")"),
         adj = c(0, 0.5), cex = 0.9)
    y_pos <- y_pos - item_height
  }
}
dev.off()
message("Saved cross-lagged network plots: net_with_self_loops.png, net_without_self_loops.png")

# 9. CENTRALITY PLOTS (In/Out Expected Influence)
plot_in <- centralityPlot(
  list("Week 0" = nwc1),
  weighted = TRUE,
  signed = TRUE,
  scale = "z-scores",
  labels = colnames_abb,
  include = "InExpectedInfluence",
  orderBy = "InExpectedInfluence"
)

plot_out <- centralityPlot(
  list("Week 0" = nwc1),
  weighted = TRUE,
  signed = TRUE,
  scale = "z-scores",
  labels = colnames_abb,
  include = "OutExpectedInfluence",
  orderBy = "OutExpectedInfluence"
)

combined_plot <- plot_in$plot + plot_out$plot + plot_layout(ncol = 2)
ggsave(
  filename = file.path(save_path_cross_lagged, "cross_lagged_centrality.png"),
  plot = combined_plot,
  width = 10, height = 10, units = "in", dpi = 300
)

# 10. BOOTSTRAP ANALYSIS
# Prepare data for bootnet (no covariates in bootstrap data)
boot_data <- data.frame(w0b, w8b, ipw = weights_global)

# Custom estimator for cross-lagged network
CLPN.OR <- function(data) {
  k <- 10
  ipw <- data[, "ipw"]
  adj <- matrix(0, k, k)
  for (i in 1:k) {
    fit <- cv.glmnet(
      x           = as.matrix(data[, 1:k]),
      y           = data[, k + i],
      family      = "gaussian",
      alpha       = 1,
      standardize = TRUE,
      weights     = ipw
    )
    lam <- fit$lambda.min
    adj[, i] <- coef(fit, s = lam)[-1]  # skip intercept
  }
  return(adj)
}

# Estimate network object
net.2 <- estimateNetwork(boot_data, fun = CLPN.OR, labels = colnames_abb, directed = TRUE)

# Plot estimated network
png(file.path(save_path_cross_lagged, "net_average.png"),
    width = 14, height = 8, units = "in", res = 300)
par(fig = c(0, 0.65, 0, 1), mar = c(2, 4, 2, 2))
plot(net.2, groups = groups_vec, legend = FALSE, color = node_colors)

par(fig = c(0.65, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0))
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 3), ylim = c(0, 6))
x_left <- 0; x_right <- 2.5; y_pos <- 5.8; item_height <- 0.2
for (group_name in names(groups)) {
  text(x_left + 0.05, y_pos - 0.18, group_name, font = 2, cex = 1.5, adj = c(0, 0.5))
  y_pos <- y_pos - item_height * 2
  for (item in groups[[group_name]]) {
    rect(x_left, y_pos - item_height, x_right, y_pos,
         col = group_colors[group_name], border = "black")
    text(x_left + 0.05, y_pos - item_height / 2,
         paste0(item, " (", item_shortnames[item], ")"),
         adj = c(0, 0.5), cex = 0.9)
    y_pos <- y_pos - item_height
  }
}
dev.off()

# 11. RUN BOOTSTRAP
set.seed(123)
nonParBoot.w2 <- bootnet(
  net.2,
  type = "nonparametric",
  nBoots = 1000,
  directed = TRUE,
  statistics = c("edge", "outExpectedInfluence", "inExpectedInfluence"),
  nCores = 4
)

set.seed(123)
caseBoot.w2 <- bootnet(
  net.2,
  type = "case",
  nBoots = 1000,
  directed = TRUE,
  statistics = c("outExpectedInfluence", "inExpectedInfluence")
)

# 12. EXTRACT AND SAVE RESULTS
cs_results <- corStability(caseBoot.w2)
boot_results <- nonParBoot.w2$sampleTable
boot_table <- nonParBoot.w2$bootTable

# Compute bootstrap statistics per edge
boot_stats <- boot_table %>%
  mutate(edge_name = paste(node1, "->", node2)) %>%
  group_by(edge_name) %>%
  summarise(
    se = sd(value, na.rm = TRUE),
    ci_lower = quantile(value, 0.025, na.rm = TRUE),
    ci_upper = quantile(value, 0.975, na.rm = TRUE),
    p_value = 2 * pmin(mean(value > 0, na.rm = TRUE), 1 - mean(value > 0, na.rm = TRUE)),
    .groups = 'drop'
  )

# Combine with sample estimates
results_table <- boot_results %>%
  mutate(edge_name = paste(node1, "->", node2)) %>%
  select(from = node1, to = node2, beta = value, type, name, edge_name) %>%
  left_join(boot_stats, by = "edge_name") %>%
  mutate(
    cs_inEI = cs_results["inExpectedInfluence"],
    cs_outEI = cs_results["outExpectedInfluence"]
  ) %>%
  select(-edge_name)

# Save to Excel
write.xlsx(results_table, file = file.path(save_path_cross_lagged, "network_path_coefficients_with_stability.xlsx"))

# 13. SAVE ADDITIONAL PLOTS
ggsave(
  filename = file.path(save_path_cross_lagged, "nonParBoot_w2_sample.png"),
  plot = plot(nonParBoot.w2, order = "sample", labels = FALSE),
  width = 6, height = 8, dpi = 300
)

ggsave(
  filename = file.path(save_path_cross_lagged, "caseBoot_w2_centrality.png"),
  plot = plot(caseBoot.w2, statistics = c("outExpectedInfluence", "inExpectedInfluence")),
  width = 8, height = 8, dpi = 300
)

ggsave(
  filename = file.path(save_path_cross_lagged, "nonParBoot_w2_outEI_diff.png"),
  plot = plot(nonParBoot.w2, "outExpectedInfluence", plot = "difference", order = "sample"),
  width = 8, height = 8, dpi = 300
)

ggsave(
  filename = file.path(save_path_cross_lagged, "nonParBoot_w2_inEI_diff.png"),
  plot = plot(nonParBoot.w2, "inExpectedInfluence", plot = "difference", order = "sample"),
  width = 8, height = 8, dpi = 300
)

ggsave(
  filename = file.path(save_path_cross_lagged, "nonParBoot_w2_edge_diff.png"),
  plot = plot(nonParBoot.w2, "edge", plot = "difference", onlyNonZero = TRUE, order = "sample"),
  width = 8, height = 8, dpi = 300
)

# 14. EXPORT TOP EDGES
diag(adjMat4) <- 1  # temporarily restore for ranking
res.w2 <- order(adjMat4, decreasing = TRUE)[seq_len(400)]
pos.w2 <- arrayInd(res.w2, dim(adjMat4), useNames = TRUE)
posWithLabs.w2 <- data.frame(
  nodeOut   = colnames_abb[pos.w2[, 1]],
  nodeIn    = colnames_abb[pos.w2[, 2]],
  value     = adjMat4[res.w2],
  DomainOut = groups_vec[pos.w2[, 1]],
  DomainIn  = groups_vec[pos.w2[, 2]]
)
edges <- data.frame(posWithLabs.w2)

# Save edge list
write_xlsx(edges, path = file.path("results", "edges_new.xlsx"))
write_xlsx(
  list("Global" = edges),
  path = file.path("results", "edges_cross_lagged_subgroups.xlsx")
)

message("Cross-lagged analysis completed. All results saved to:", save_path_cross_lagged)