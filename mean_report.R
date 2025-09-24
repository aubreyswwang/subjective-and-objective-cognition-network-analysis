## Cross-sectional and longitudinal network analysis

# 1. LOAD REQUIRED PACKAGES
pkgs <- c("readxl", "dplyr", "haven", "psych", "car", "NetworkComparisonTest",
          "ggplot2", "mgm", "mice", "bootnet", "qgraph", "psychTools", 
          "glmnet", "lavaan", "readr", "EstimateGroupNetwork", "networktools", 
          "NetworkToolbox", "Hmisc", "wCorr", "writexl", "patchwork", "openxlsx",
          "huge", "officer", "flextable", "gridExtra", "tidyr", "forcats", "Rgraphviz")

# Install and load packages if not already available
invisible(lapply(pkgs, function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}))

# 2. READ DATA FILE
network_data <- read_excel(
  "mean_imputed_with_single_ipw.xlsx",
  sheet = "Mean_Imputed"
)

# Keep only retained participants
network_data <- network_data[network_data$Retained == "retained" & !is.na(network_data$Retained), ]

# Extract pre- and post-treatment variables (columns 13–22 and 36–45)
dataqids_biomarkers_pre  <- network_data[, 13:22]
dataqids_biomarkers_post <- network_data[, 36:45]

# Subgroup data: subgroup == 0 (remission), subgroup == 1 (non-remission)
subgroup_rem_pre      <- network_data[network_data$subgroup == 0 & !is.na(network_data$subgroup), 13:22]
subgroup_nonrem_pre   <- network_data[network_data$subgroup == 1 & !is.na(network_data$subgroup), 13:22]
subgroup_rem_post     <- network_data[network_data$subgroup == 0 & !is.na(network_data$subgroup), 36:45]
subgroup_nonrem_post  <- network_data[network_data$subgroup == 1 & !is.na(network_data$subgroup), 36:45]

# Extract inverse probability weights (IPW)
weights_global   <- network_data$IPW_w[network_data$Retained == "retained" & !is.na(network_data$Retained)]
weights_rem      <- network_data$IPW_w[network_data$subgroup == 0 & !is.na(network_data$subgroup)]
weights_nonrem   <- network_data$IPW_w[network_data$subgroup == 1 & !is.na(network_data$subgroup)]

# Define output directories (relative paths)
save_path           <- "results"
save_path_net       <- file.path(save_path, "NET")
save_path_centrality<- file.path(save_path, "centrality")
save_path_stability <- file.path(save_path, "stability")
save_path_accuracy  <- file.path(save_path, "accuracy")
save_path_diff_test <- file.path(save_path, "difference_test")

# Create directories if they don't exist
dir.create(save_path,           showWarnings = FALSE, recursive = TRUE)
dir.create(save_path_net,       showWarnings = FALSE)
dir.create(save_path_centrality,showWarnings = FALSE)
dir.create(save_path_stability, showWarnings = FALSE)
dir.create(save_path_accuracy,  showWarnings = FALSE)
dir.create(save_path_diff_test, showWarnings = FALSE)

# 3. DEFINE VARIABLE NAMES AND GROUPS
colnames_abb <- c(
  "S-AT", "S-ME", "S-EF",
  "O-ME", "O-EF", "O-AT", "O-PS",
  "SM", "GI", "EN"
)

items_clean <- c(
  "Subjective Attention", "Subjective Memory", "Subjective Executive Functions",
  "Objective Memory", "Objective Executive Functions", "Objective Attention", "Objective Processing Speed",
  "Sad Mood", "General Interest", "Energy"
)

groups <- list(
  "Subjective Cognition" = c(1:3),
  "Objective Cognition"  = c(4:7),
  "Depressive Symptom"   = c(8:10)
)

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
item_shortnames <- setNames(colnames_abb, items_clean)
groups_named <- split(items_clean, groups_vec)

# Apply consistent column names to all datasets
data_list_pre <- list(
  dataqids_biomarkers_pre = dataqids_biomarkers_pre,
  subgroup_rem_pre = subgroup_rem_pre,
  subgroup_nonrem_pre = subgroup_nonrem_pre
)

data_list_post <- list(
  dataqids_biomarkers_post = dataqids_biomarkers_post,
  subgroup_rem_post = subgroup_rem_post,
  subgroup_nonrem_post = subgroup_nonrem_post
)

data_list_pre  <- lapply(data_list_pre,  function(df) { colnames(df) <- colnames_abb; df })
data_list_post <- lapply(data_list_post, function(df) { colnames(df) <- colnames_abb; df })

# Export to global environment (as in original script)
list2env(data_list_pre,  envir = .GlobalEnv)
list2env(data_list_post, envir = .GlobalEnv)

# 4. SELECT VARIABLES AND CLEAN DATA
biomarkers_qids_pre      <- as.data.frame(dataqids_biomarkers_pre)[, 1:10, drop = FALSE]
biomarkers_qids_post     <- as.data.frame(dataqids_biomarkers_post)[, 1:10, drop = FALSE]
biomarkers_rem_pre       <- as.data.frame(subgroup_rem_pre)[, 1:10, drop = FALSE]
biomarkers_nonrem_pre    <- as.data.frame(subgroup_nonrem_pre)[, 1:10, drop = FALSE]
biomarkers_rem_post      <- as.data.frame(subgroup_rem_post)[, 1:10, drop = FALSE]
biomarkers_nonrem_post   <- as.data.frame(subgroup_nonrem_post)[, 1:10, drop = FALSE]

# Convert to numeric and remove rows with missing values
clean_data <- function(df) {
  df %>%
    mutate(across(everything(), ~ as.numeric(as.character(.)))) %>%
    na.omit()
}

biomarkers_qids_pre      <- clean_data(biomarkers_qids_pre)
biomarkers_qids_post     <- clean_data(biomarkers_qids_post)
biomarkers_rem_pre       <- clean_data(biomarkers_rem_pre)
biomarkers_nonrem_pre    <- clean_data(biomarkers_nonrem_pre)
biomarkers_rem_post      <- clean_data(biomarkers_rem_post)
biomarkers_nonrem_post   <- clean_data(biomarkers_nonrem_post)

# Nonparanormal transformation for non-normal data
biomarkers_qids_pre      <- huge.npn(biomarkers_qids_pre)
biomarkers_qids_post     <- huge.npn(biomarkers_qids_post)
biomarkers_rem_pre       <- huge.npn(biomarkers_rem_pre)
biomarkers_nonrem_pre    <- huge.npn(biomarkers_nonrem_pre)
biomarkers_rem_post      <- huge.npn(biomarkers_rem_post)
biomarkers_nonrem_post   <- huge.npn(biomarkers_nonrem_post)

# 5. ESTIMATE MGM NETWORKS
p1 <- ncol(biomarkers_qids_pre)
nodepredict_biomarkers_qids_pre <- mgm(
  data = biomarkers_qids_pre,
  type = rep('g', p1),
  level = rep(1, p1),
  lambdaSel = "CV",
  ruleReg = "OR"
)
predictabil_biomarkers_qids_pre <- predict(
  object = nodepredict_biomarkers_qids_pre,
  data = biomarkers_qids_pre,
  errorCon = "R2"
)

p2 <- ncol(biomarkers_qids_post)
nodepredict_biomarkers_qids_post <- mgm(
  data = biomarkers_qids_post,
  type = rep('g', p2),
  level = rep(1, p2),
  lambdaSel = "CV",
  ruleReg = "OR"
)
predictabil_biomarkers_qids_post <- predict(
  object = nodepredict_biomarkers_qids_post,
  data = biomarkers_qids_post,
  errorCon = "R2"
)

p3 <- ncol(biomarkers_rem_pre)
nodepredict_biomarkers_rem_pre <- mgm(
  data = biomarkers_rem_pre,
  type = rep('g', p3),
  level = rep(1, p3),
  lambdaSel = "CV",
  ruleReg = "OR"
)
predictabil_biomarkers_rem_pre <- predict(
  object = nodepredict_biomarkers_rem_pre,
  data = biomarkers_rem_pre,
  errorCon = "R2"
)

p4 <- ncol(biomarkers_nonrem_pre)
nodepredict_biomarkers_nonrem_pre <- mgm(
  data = biomarkers_nonrem_pre,
  type = rep('g', p4),
  level = rep(1, p4),
  lambdaSel = "CV",
  ruleReg = "OR"
)
predictabil_biomarkers_nonrem_pre <- predict(
  object = nodepredict_biomarkers_nonrem_pre,
  data = biomarkers_nonrem_pre,
  errorCon = "R2"
)

p5 <- ncol(biomarkers_rem_post)
nodepredict_biomarkers_rem_post <- mgm(
  data = biomarkers_rem_post,
  type = rep('g', p5),
  level = rep(1, p5),
  lambdaSel = "CV",
  ruleReg = "OR"
)
predictabil_biomarkers_rem_post <- predict(
  object = nodepredict_biomarkers_rem_post,
  data = biomarkers_rem_post,
  errorCon = "R2"
)

p6 <- ncol(biomarkers_nonrem_post)
nodepredict_biomarkers_nonrem_post <- mgm(
  data = biomarkers_nonrem_post,
  type = rep('g', p6),
  level = rep(1, p6),
  lambdaSel = "CV",
  ruleReg = "OR"
)
predictabil_biomarkers_nonrem_post <- predict(
  object = nodepredict_biomarkers_nonrem_post,
  data = biomarkers_nonrem_post,
  errorCon = "R2"
)

# 5.1 Calculate mean predictability (R²)
mean_predict_biomarkers_qids_pre      <- round(mean(predictabil_biomarkers_qids_pre$error$R2), digits = 2)
mean_predict_biomarkers_qids_post     <- round(mean(predictabil_biomarkers_qids_post$error$R2), digits = 2)
mean_predict_biomarkers_rem_pre       <- round(mean(predictabil_biomarkers_rem_pre$error$R2), digits = 2)
mean_predict_biomarkers_nonrem_pre    <- round(mean(predictabil_biomarkers_nonrem_pre$error$R2), digits = 2)
mean_predict_biomarkers_rem_post      <- round(mean(predictabil_biomarkers_rem_post$error$R2), digits = 2)
mean_predict_biomarkers_nonrem_post   <- round(mean(predictabil_biomarkers_nonrem_post$error$R2), digits = 2)

# Save R² results
results <- c(
  "mean_predict_biomarkers_qids_pre"      = mean_predict_biomarkers_qids_pre,
  "mean_predict_biomarkers_qids_post"     = mean_predict_biomarkers_qids_post,
  "mean_predict_biomarkers_rem_pre"       = mean_predict_biomarkers_rem_pre,
  "mean_predict_biomarkers_nonrem_pre"    = mean_predict_biomarkers_nonrem_pre,
  "mean_predict_biomarkers_rem_post"      = mean_predict_biomarkers_rem_post,
  "mean_predict_biomarkers_nonrem_post"   = mean_predict_biomarkers_nonrem_post
)

writeLines("Mean R2 Values:\n", con = file.path(save_path_centrality, "R2_results.txt"))
write.table(results, file = file.path(save_path_centrality, "R2_results.txt"),
            append = TRUE, sep = " = ", col.names = FALSE, row.names = TRUE, quote = FALSE)
cat("R2 results saved to:", save_path_centrality, "\n")

# 7. CALCULATE WEIGHTED EBICGLASSO NETWORKS
cov_pre  <- cov.wt(biomarkers_qids_pre,  wt = weights_global,  cor = FALSE)$cov
cov_post <- cov.wt(biomarkers_qids_post, wt = weights_global, cor = FALSE)$cov
cov_rem_pre  <- cov.wt(biomarkers_rem_pre,  wt = weights_rem,  cor = FALSE)$cov
cov_nonrem_pre  <- cov.wt(biomarkers_nonrem_pre,  wt = weights_nonrem,  cor = FALSE)$cov
cov_rem_post <- cov.wt(biomarkers_rem_post, wt = weights_rem, cor = FALSE)$cov
cov_nonrem_post <- cov.wt(biomarkers_nonrem_post, wt = weights_nonrem, cor = FALSE)$cov

cormatrix_weighted_pre  <- cov2cor(cov_pre)
cormatrix_weighted_post <- cov2cor(cov_post)
cormatrix_weighted_rem_pre  <- cov2cor(cov_rem_pre)
cormatrix_weighted_nonrem_pre  <- cov2cor(cov_nonrem_pre)
cormatrix_weighted_rem_post <- cov2cor(cov_rem_post)
cormatrix_weighted_nonrem_post <- cov2cor(cov_nonrem_post)

# Effective sample size (ESS)
ess_global  <- sum(weights_global)^2  / sum(weights_global^2)
ess_rem     <- sum(weights_rem)^2     / sum(weights_rem^2)
ess_nonrem  <- sum(weights_nonrem)^2  / sum(weights_nonrem^2)

# Network plotting function
plot_qgraph_with_legend <- function(cor_matrix, predictability, labels,
                                    groups_vec, groups, group_colors,
                                    node_colors, sample_size,
                                    filename, title = "") {
  png(filename = file.path(save_path_net, filename), width = 8, height = 8.5, units = "in", res = 300)
  par(mar = c(2, 2, 3, 2), xpd = NA)
  net_obj <- qgraph(
    cor_matrix,
    graph = "glasso",
    gamma = 0.5,
    layout = "spring",
    groups = groups_vec,
    sampleSize = sample_size,
    labels = labels,
    color = node_colors,
    pie = as.numeric(as.character(predictability$errors[, 2])),
    pieColor = rep('#377EB8', length(labels)),
    border.color = "black",
    border.width = 1,
    theme = "Borkulo",
    legend = FALSE,
    vsize = 8.5,
    edge.labels = TRUE,
    edge.label.cex = 0.7,
    fade = TRUE,
    title = NULL
  )
  title(main = title, col.main = "black", font.main = 2, cex.main = 2.0, adj = 0, line = 0.3)
  dev.off()
  return(net_obj)
}

# Generate all network plots
net_0 <- plot_qgraph_with_legend(
  cor_matrix = cormatrix_weighted_pre,
  predictability = predictabil_biomarkers_qids_pre,
  labels = colnames_abb,
  groups_vec = groups_vec,
  groups = groups_named,
  group_colors = group_colors,
  node_colors = node_colors,
  sample_size = ess_global,
  filename = "net_0.png",
  title = "Week 0"
)

net_8 <- plot_qgraph_with_legend(
  cor_matrix = cormatrix_weighted_post,
  predictability = predictabil_biomarkers_qids_post,
  labels = colnames_abb,
  groups_vec = groups_vec,
  groups = groups_named,
  group_colors = group_colors,
  node_colors = node_colors,
  sample_size = ess_global,
  filename = "net_8.png",
  title = "Week 8"
)

net_rem_0 <- plot_qgraph_with_legend(
  cor_matrix = cormatrix_weighted_rem_pre,
  predictability = predictabil_biomarkers_rem_pre,
  labels = colnames_abb,
  groups_vec = groups_vec,
  groups = groups_named,
  group_colors = group_colors,
  node_colors = node_colors,
  sample_size = ess_rem,
  filename = "net_rem_0.png",
  title = "Week 0_Remission"
)

net_nonrem_0 <- plot_qgraph_with_legend(
  cor_matrix = cormatrix_weighted_nonrem_pre,
  predictability = predictabil_biomarkers_nonrem_pre,
  labels = colnames_abb,
  groups_vec = groups_vec,
  groups = groups_named,
  group_colors = group_colors,
  node_colors = node_colors,
  sample_size = ess_nonrem,
  filename = "net_nonrem_0.png",
  title = "Week 0_Non-remission"
)

net_rem_8 <- plot_qgraph_with_legend(
  cor_matrix = cormatrix_weighted_rem_post,
  predictability = predictabil_biomarkers_rem_post,
  labels = colnames_abb,
  groups_vec = groups_vec,
  groups = groups_named,
  group_colors = group_colors,
  node_colors = node_colors,
  sample_size = ess_rem,
  filename = "net_rem_8.png",
  title = "Week 8_Remission"
)

net_nonrem_8 <- plot_qgraph_with_legend(
  cor_matrix = cormatrix_weighted_nonrem_post,
  predictability = predictabil_biomarkers_nonrem_post,
  labels = colnames_abb,
  groups_vec = groups_vec,
  groups = groups_named,
  group_colors = group_colors,
  node_colors = node_colors,
  sample_size = ess_nonrem,
  filename = "net_nonrem_8.png",
  title = "Week 8_Non-remission"
)

# Save custom legends
save_single_group_legend <- function(group_name, items, group_color, item_shortnames, filename, width = 3.0, res = 120) {
  n_items <- length(items)
  item_height     <- 0.3
  title_gap       <- 0.25
  bottom_padding  <- 0.10
  top_padding     <- 0.10
  text_y_offset   <- 0.01
  group_name_cex  <- 1.5
  item_text_cex   <- 1.2
  bar_width       <- 2.9
  total_height <- top_padding + 0.2 + title_gap + n_items * item_height + bottom_padding
  total_width  <- width
  png(file.path(save_path_net, filename), total_width, total_height, "in", res = res)
  par(mar = c(0.03, 0.03, 0, 0.03), xpd = NA)
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 3), ylim = c(0, total_height))
  y_pos <- total_height - top_padding - 0.2
  text(0.05, y_pos, group_name, font = 2, cex = group_name_cex, adj = c(0, 0.5))
  y_pos <- y_pos - title_gap
  for (item in items) {
    rect(0.015, y_pos - item_height, 0.05 + bar_width, y_pos, col = group_color, border = "black")
    text(0.05, y_pos - item_height / 2 + text_y_offset,
         paste0(item, " (", item_shortnames[item], ")"),
         adj = c(0, 0.5), cex = item_text_cex)
    y_pos <- y_pos - item_height
  }
  dev.off()
}

for (grp_name in names(groups_named)) {
  save_single_group_legend(
    group_name = grp_name,
    items = groups_named[[grp_name]],
    group_color = group_colors[grp_name],
    item_shortnames = item_shortnames,
    filename = paste0("legend_", gsub(" ", "_", grp_name), ".png"),
    width = 4,
    res = 120
  )
}

# Extract and save edge matrices
extract_edge_matrix <- function(qgraph_obj, network_name) {
  weight_matrix <- getWmat(qgraph_obj)
  df <- as.data.frame(weight_matrix)
  rownames(df) <- colnames(df) <- colnames_abb
  df$Source_Network <- network_name
  return(df)
}

all_edge_matrices <- list()
all_edge_matrices[["net_0"]]               <- extract_edge_matrix(net_0, "Week 0")
all_edge_matrices[["net_8"]]               <- extract_edge_matrix(net_8, "Week 8")
all_edge_matrices[["net_rem_0"]]           <- extract_edge_matrix(net_rem_0, "Week 0_Remission")
all_edge_matrices[["net_nonrem_0"]]        <- extract_edge_matrix(net_nonrem_0, "Week 0_Non-remission")
all_edge_matrices[["net_rem_8"]]           <- extract_edge_matrix(net_rem_8, "Week 8_Remission")
all_edge_matrices[["net_nonrem_8"]]        <- extract_edge_matrix(net_nonrem_8, "Week 8_Non-remission")

wb_edge <- createWorkbook()
for (sheet_name in names(all_edge_matrices)) {
  addWorksheet(wb_edge, sheet_name)
  writeData(wb_edge, sheet = sheet_name, x = all_edge_matrices[[sheet_name]])
}
saveWorkbook(wb_edge, file = file.path(save_path_net, "All_Network_Edge_Matrices.xlsx"), overwrite = TRUE)

# 8. CENTRALITY ESTIMATES
save_centrality_plot <- function(filename, graph_list, labels, include, orderBy) {
  png(filename = file.path(save_path_centrality, filename), width = 6, height = 8, units = "in", res = 300)
  result <- centralityPlot(
    graph_list,
    weighted = TRUE,
    signed = TRUE,
    scale = "z-scores",
    labels = labels,
    include = include,
    orderBy = orderBy
  )
  print(result$plot)
  dev.off()
  return(result)
}

# Expected Influence
centrality_biomarkers_qids_EI <- save_centrality_plot(
  "centrality_EI_0_vs_8.png",
  list("Week 0" = net_0, "Week 8" = net_8),
  labels = colnames_abb,
  include = "ExpectedInfluence",
  orderBy = "ExpectedInfluence"
)

centrality_biomarkers_qids_EI_rem <- save_centrality_plot(
  "centrality_EI_rem.png",
  list("Week 0_Remission" = net_rem_0, "Week 8_Remission" = net_rem_8),
  labels = colnames_abb,
  include = "ExpectedInfluence",
  orderBy = "ExpectedInfluence"
)

centrality_biomarkers_qids_EI_nonrem <- save_centrality_plot(
  "centrality_EI_nonrem.png",
  list("Week 0_Non-remission" = net_nonrem_0, "Week 8_Non-remission" = net_nonrem_8),
  labels = colnames_abb,
  include = "ExpectedInfluence",
  orderBy = "ExpectedInfluence"
)

# Strength
centrality_biomarkers_qids_strength <- save_centrality_plot(
  "centrality_Strength_0_vs_8.png",
  list("Week 0" = net_0, "Week 8" = net_8),
  labels = colnames_abb,
  include = "Strength",
  orderBy = "Strength"
)

centrality_biomarkers_qids_strength_rem <- save_centrality_plot(
  "centrality_Strength_rem.png",
  list("Week 0_Remission" = net_rem_0, "Week 8_Remission" = net_rem_8),
  labels = colnames_abb,
  include = "Strength",
  orderBy = "Strength"
)

centrality_biomarkers_qids_strength_nonrem <- save_centrality_plot(
  "centrality_Strength_nonrem.png",
  list("Week 0_Non-remission" = net_nonrem_0, "Week 8_Non-remission" = net_nonrem_8),
  labels = colnames_abb,
  include = "Strength",
  orderBy = "Strength"
)

# Combined
centrality_combined <- save_centrality_plot(
  "centrality_combined_0_vs_8.png",
  list("Week 0" = net_0, "Week 8" = net_8),
  labels = colnames_abb,
  include = c("Strength", "ExpectedInfluence"),
  orderBy = "ExpectedInfluence"
)

centrality_combined_rem <- save_centrality_plot(
  "centrality_combined_rem.png",
  list("Week 0_Remission" = net_rem_0, "Week 8_Remission" = net_rem_8),
  labels = colnames_abb,
  include = c("Strength", "ExpectedInfluence"),
  orderBy = "ExpectedInfluence"
)

centrality_combined_nonrem <- save_centrality_plot(
  "centrality_combined_nonrem.png",
  list("Week 0_Non-remission" = net_nonrem_0, "Week 8_Non-remission" = net_nonrem_8),
  labels = colnames_abb,
  include = c("Strength", "ExpectedInfluence"),
  orderBy = "ExpectedInfluence"
)

# Save centrality results to Excel
all_results <- list()
all_results[["EI_w0vs8"]]            <- cbind(Analysis = "EI_w0vs8", centrality_biomarkers_qids_EI$data)
all_results[["EI_w0vs8_Rem"]]        <- cbind(Analysis = "EI_w0vs8_Rem", centrality_biomarkers_qids_EI_rem$data)
all_results[["EI_w0vs8_NonRem"]]     <- cbind(Analysis = "EI_w0vs8_NonRem", centrality_biomarkers_qids_EI_nonrem$data)
all_results[["Strength_w0vs8"]]      <- cbind(Analysis = "Strength_w0vs8", centrality_biomarkers_qids_strength$data)
all_results[["Strength_w0vs8_Rem"]]  <- cbind(Analysis = "Strength_w0vs8_Rem", centrality_biomarkers_qids_strength_rem$data)
all_results[["Strength_w0vs8_NonRem"]] <- cbind(Analysis = "Strength_w0vs8_NonRem", centrality_biomarkers_qids_strength_nonrem$data)

combined_all <- bind_rows(all_results)
wb2 <- createWorkbook()
addWorksheet(wb2, "Centrality_Results")
writeData(wb2, "Centrality_Results", combined_all)
saveWorkbook(wb2, file = file.path(save_path_centrality, "Combined_Centrality_Results.xlsx"), overwrite = TRUE)

# 10. NETWORK COMPARISON TESTS (NCT)
myWeightedGlasso <- function(data, weights, gamma, verbose = FALSE, ...) {
  covm <- cov.wt(data, wt = weights, cor = FALSE)$cov
  C    <- cov2cor(covm)
  net  <- qgraph::EBICglasso(C, n = length(weights), gamma = gamma)
  return(net)
}

set.seed(2025)
nct_biomarkers <- NCT(
  data1         = as.data.frame(biomarkers_qids_pre),
  data2         = as.data.frame(biomarkers_qids_post),
  estimator     = myWeightedGlasso,
  estimatorArgs = list(weights = weights_global, gamma = 0.5),
  it               = 1000,
  test.edges       = TRUE,
  edges            = "all",
  p.adjust.methods = "holm",
  test.centrality  = TRUE,
  centrality       = "expectedInfluence"
)

set.seed(2025)
nct_biomarkers_rem <- NCT(
  data1         = as.data.frame(biomarkers_rem_pre),
  data2         = as.data.frame(biomarkers_rem_post),
  estimator     = myWeightedGlasso,
  estimatorArgs = list(weights = weights_rem, gamma = 0.5),
  it               = 1000,
  test.edges       = TRUE,
  edges            = "all",
  p.adjust.methods = "holm",
  test.centrality  = TRUE,
  centrality       = "expectedInfluence"
)

set.seed(2025)
nct_biomarkers_nonrem <- NCT(
  data1         = as.data.frame(biomarkers_nonrem_pre),
  data2         = as.data.frame(biomarkers_nonrem_post),
  estimator     = myWeightedGlasso,
  estimatorArgs = list(weights = weights_nonrem, gamma = 0.5),
  it               = 1000,
  test.edges       = TRUE,
  edges            = "all",
  p.adjust.methods = "holm",
  test.centrality  = TRUE,
  centrality       = "expectedInfluence"
)

set.seed(2025)
nct_rem_vs_nonrem_pre <- NCT(
  data1         = as.data.frame(biomarkers_rem_pre),
  data2         = as.data.frame(biomarkers_nonrem_pre),
  estimator = function(data, ...) {
    if (nrow(data) == nrow(biomarkers_rem_pre) && all(colnames(data) == colnames(biomarkers_rem_pre))) {
      weights <- weights_rem
    } else if (nrow(data) == nrow(biomarkers_nonrem_pre) && all(colnames(data) == colnames(biomarkers_nonrem_pre))) {
      weights <- weights_nonrem
    } else {
      stop("Unknown group")
    }
    myWeightedGlasso(data, weights = weights, gamma = 0.5)
  },
  estimatorArgs = list(),
  it            = 1000,
  test.edges    = TRUE,
  edges         = "all",
  test.centrality = TRUE,
  centrality      = "expectedInfluence",
  p.adjust.methods = "holm"
)

set.seed(2025)
nct_rem_vs_nonrem_post <- NCT(
  data1         = as.data.frame(biomarkers_rem_post),
  data2         = as.data.frame(biomarkers_nonrem_post),
  estimator = function(data, ...) {
    if (nrow(data) == nrow(biomarkers_rem_post) && all(colnames(data) == colnames(biomarkers_rem_post))) {
      weights <- weights_rem
    } else if (nrow(data) == nrow(biomarkers_nonrem_post) && all(colnames(data) == colnames(biomarkers_nonrem_post))) {
      weights <- weights_nonrem
    } else {
      stop("Unknown group")
    }
    myWeightedGlasso(data, weights = weights, gamma = 0.5)
  },
  estimatorArgs = list(),
  it            = 1000,
  test.edges    = TRUE,
  edges         = "all",
  test.centrality = TRUE,
  centrality      = "expectedInfluence",
  p.adjust.methods = "holm"
)

# Export NCT results to Excel
to_text <- function(x) {
  if (is.null(x)) return("NULL")
  paste(capture.output(print(x)), collapse = "\n")
}

create_summary_df <- function(nct_obj) {
  data.frame(
    Metric = c(
      "Global strength difference",
      "Global strength per network",
      "p-value (global strength)",
      "Maximum edge difference",
      "p-value (max edge)"
    ),
    Value = sapply(list(
      nct_obj$glstrinv.real,
      nct_obj$glstrinv.sep,
      nct_obj$glstrinv.pval,
      nct_obj$nwinv.real,
      nct_obj$nwinv.pval
    ), to_text),
    stringsAsFactors = FALSE
  )
}

export_nct_to_excel <- function(save_path, file_name = "nct_summary_output_raw.xlsx") {
  wb <- createWorkbook()
  add_result_sheet <- function(sheet_name, nct_obj) {
    addWorksheet(wb, sheet_name)
    row_pos <- 1
    summary_df <- create_summary_df(nct_obj)
    writeData(wb, sheet = sheet_name, x = summary_df, startRow = row_pos)
    row_pos <- row_pos + nrow(summary_df) + 2
    writeData(wb, sheet = sheet_name, x = "Edge difference matrix (einv.real)", startRow = row_pos)
    row_pos <- row_pos + 1
    writeData(wb, sheet = sheet_name, x = nct_obj$einv.real, startRow = row_pos, rowNames = TRUE)
    row_pos <- row_pos + nrow(nct_obj$einv.real) + 2
    writeData(wb, sheet = sheet_name, x = "Edge p-values (einv.pvals)", startRow = row_pos)
    row_pos <- row_pos + 1
    writeData(wb, sheet = sheet_name, x = nct_obj$einv.pvals, startRow = row_pos, rowNames = TRUE)
    row_pos <- row_pos + nrow(nct_obj$einv.pvals) + 2
    writeData(wb, sheet = sheet_name, x = "Centrality differences (diffcen.real)", startRow = row_pos)
    row_pos <- row_pos + 1
    writeData(wb, sheet = sheet_name,
              x = data.frame(Node = row.names(nct_obj$diffcen.real), Difference = nct_obj$diffcen.real),
              startRow = row_pos)
    row_pos <- row_pos + length(nct_obj$diffcen.real) + 2
    writeData(wb, sheet = sheet_name, x = "Centrality p-values (diffcen.pval)", startRow = row_pos)
    row_pos <- row_pos + 1
    writeData(wb, sheet = sheet_name,
              x = data.frame(Node = row.names(nct_obj$diffcen.pval), p_value = nct_obj$diffcen.pval),
              startRow = row_pos)
  }
  add_result_sheet("Global_pre_vs_post", nct_biomarkers)
  add_result_sheet("REM_pre_vs_post", nct_biomarkers_rem)
  add_result_sheet("NonREM_pre_vs_post", nct_biomarkers_nonrem)
  add_result_sheet("REM_vs_NonREM_at_pre", nct_rem_vs_nonrem_pre)
  add_result_sheet("REM_vs_NonREM_at_post", nct_rem_vs_nonrem_post)
  saveWorkbook(wb, file = file.path(save_path, file_name), overwrite = TRUE)
  cat("NCT Excel file saved to:", file.path(save_path, file_name), "\n")
}
export_nct_to_excel(save_path)

# 11. ACCURACY & STABILITY ANALYSES
complete_cases <- complete.cases(biomarkers_qids_pre, biomarkers_qids_post, weights_global)
complete_cases_rem <- complete.cases(biomarkers_rem_pre, biomarkers_rem_post, weights_rem)
complete_cases_nonrem <- complete.cases(biomarkers_nonrem_pre, biomarkers_nonrem_post, weights_nonrem)

df_pre <- as.matrix(biomarkers_qids_pre[complete_cases, ])
df_post <- as.matrix(biomarkers_qids_post[complete_cases, ])
weights <- as.numeric(weights_global[complete_cases])
df_pre_rem <- as.matrix(biomarkers_rem_pre[complete_cases_rem, ])
df_pre_nonrem <- as.matrix(biomarkers_nonrem_pre[complete_cases_nonrem, ])
df_post_rem <- as.matrix(biomarkers_rem_post[complete_cases_rem, ])
df_post_nonrem <- as.matrix(biomarkers_nonrem_post[complete_cases_nonrem, ])
weights_rem <- as.numeric(weights_rem[complete_cases_rem])
weights_nonrem <- as.numeric(weights_nonrem[complete_cases_nonrem])

# Bootstrap estimators
myWeightedGlasso_bootnet <- function(data, gamma, ...) {
  data_str <- apply(data, 1, paste, collapse = "_")
  pre_str  <- apply(df_pre, 1, paste, collapse = "_")
  post_str <- apply(df_post, 1, paste, collapse = "_")
  if (all(data_str %in% pre_str)) {
    row_ids <- match(data_str, pre_str)
  } else if (all(data_str %in% post_str)) {
    row_ids <- match(data_str, post_str)
  } else {
    stop("Data source not recognized (neither pre nor post)")
  }
  weights_local <- weights[row_ids]
  if (any(is.na(weights_local)) || length(weights_local) != nrow(data)) {
    stop("Weight matching failed!")
  }
  covm <- cov.wt(data, wt = weights_local, cor = FALSE)$cov
  C <- cov2cor(covm)
  C_pd <- as.matrix(nearPD(C, corr = TRUE)$mat)
  net <- EBICglasso(C_pd, n = length(weights_local), gamma = gamma)
  return(net)
}

myWeightedGlasso_bootnet_rem <- function(data, gamma, ...) {
  data_str <- apply(data, 1, paste, collapse = "_")
  pre_str  <- apply(df_pre_rem, 1, paste, collapse = "_")
  post_str <- apply(df_post_rem, 1, paste, collapse = "_")
  if (all(data_str %in% pre_str)) {
    row_ids <- match(data_str, pre_str)
  } else if (all(data_str %in% post_str)) {
    row_ids <- match(data_str, post_str)
  } else {
    stop("Data source not recognized (neither pre_rem nor post_rem)")
  }
  weights_local <- weights_rem[row_ids]
  if (any(is.na(weights_local)) || length(weights_local) != nrow(data)) {
    stop("Weight matching failed!")
  }
  covm <- cov.wt(data, wt = weights_local, cor = FALSE)$cov
  C <- cov2cor(covm)
  C_pd <- as.matrix(nearPD(C, corr = TRUE)$mat)
  net <- EBICglasso(C_pd, n = length(weights_local), gamma = gamma)
  return(net)
}

myWeightedGlasso_bootnet_nonrem <- function(data, gamma, ...) {
  data_str <- apply(data, 1, paste, collapse = "_")
  pre_str  <- apply(df_pre_nonrem, 1, paste, collapse = "_")
  post_str <- apply(df_post_nonrem, 1, paste, collapse = "_")
  if (all(data_str %in% pre_str)) {
    row_ids <- match(data_str, pre_str)
  } else if (all(data_str %in% post_str)) {
    row_ids <- match(data_str, post_str)
  } else {
    stop("Data source not recognized (neither pre_nonrem nor post_nonrem)")
  }
  weights_local <- weights_nonrem[row_ids]
  if (any(is.na(weights_local)) || length(weights_local) != nrow(data)) {
    stop("Weight matching failed!")
  }
  covm <- cov.wt(data, wt = weights_local, cor = FALSE)$cov
  C <- cov2cor(covm)
  C_pd <- as.matrix(nearPD(C, corr = TRUE)$mat)
  net <- EBICglasso(C_pd, n = length(weights_local), gamma = gamma)
  return(net)
}

# Case-dropping bootstrap
set.seed(1234)
casedrop_boot_pre <- bootnet(
  data = df_pre,
  default = "none",
  fun = myWeightedGlasso_bootnet,
  type = "case",
  statistics = c("edge", "expectedInfluence"),
  nBoots = 1000,
  gamma = 0.5,
  caseN = nrow(biomarkers_qids_pre),
  maxErrors = 10,
  errorHandling = "continue",
  nCores = 1
)

set.seed(1234)
casedrop_boot_post <- bootnet(
  data = df_post,
  default = "none",
  fun = myWeightedGlasso_bootnet,
  type = "case",
  statistics = c("edge", "expectedInfluence"),
  nBoots = 1000,
  gamma = 0.5,
  caseN = nrow(biomarkers_qids_post),
  maxErrors = 10,
  errorHandling = "continue",
  nCores = 1
)

set.seed(1234)
casedrop_boot_pre_rem <- bootnet(
  data = df_pre_rem,
  default = "none",
  fun = myWeightedGlasso_bootnet_rem,
  type = "case",
  statistics = c("edge", "expectedInfluence"),
  nBoots = 1000,
  gamma = 0.5,
  caseN = nrow(biomarkers_rem_pre),
  maxErrors = 10,
  errorHandling = "continue",
  nCores = 1
)

set.seed(1234)
casedrop_boot_pre_nonrem <- bootnet(
  data = df_pre_nonrem,
  default = "none",
  fun = myWeightedGlasso_bootnet_nonrem,
  type = "case",
  statistics = c("edge", "expectedInfluence"),
  nBoots = 1000,
  gamma = 0.5,
  caseN = nrow(biomarkers_nonrem_pre),
  maxErrors = 10,
  errorHandling = "continue",
  nCores = 1
)

set.seed(1234)
casedrop_boot_post_rem <- bootnet(
  data = df_post_rem,
  default = "none",
  fun = myWeightedGlasso_bootnet_rem,
  type = "case",
  statistics = c("edge", "expectedInfluence"),
  nBoots = 1000,
  gamma = 0.5,
  caseN = nrow(biomarkers_rem_post),
  maxErrors = 10,
  errorHandling = "continue",
  nCores = 1
)

set.seed(1234)
casedrop_boot_post_nonrem <- bootnet(
  data = df_post_nonrem,
  default = "none",
  fun = myWeightedGlasso_bootnet_nonrem,
  type = "case",
  statistics = c("edge", "expectedInfluence"),
  nBoots = 1000,
  gamma = 0.5,
  caseN = nrow(biomarkers_nonrem_post),
  maxErrors = 10,
  errorHandling = "continue",
  nCores = 1
)

# Save stability plots
title_theme <- theme(legend.position = "none", plot.title = element_text(size = 16, hjust = 0, face = "bold"))

ggsave(filename = file.path(save_path_stability, "EI_0w.png"),
       plot = plot(casedrop_boot_pre, statistics = "expectedInfluence", order = "sample") + title_theme + ggtitle("Expected Influence – 0w"),
       width = 8, height = 8, units = "in")

ggsave(filename = file.path(save_path_stability, "EI_8w.png"),
       plot = plot(casedrop_boot_post, statistics = "expectedInfluence", order = "sample") + title_theme + ggtitle("Expected Influence – 8w"),
       width = 8, height = 8, units = "in")

ggsave(filename = file.path(save_path_stability, "EI_0w_REM.png"),
       plot = plot(casedrop_boot_pre_rem, statistics = "expectedInfluence", order = "sample") + title_theme + ggtitle("Expected Influence – 0w (REM)"),
       width = 8, height = 8, units = "in")

ggsave(filename = file.path(save_path_stability, "EI_8w_REM.png"),
       plot = plot(casedrop_boot_post_rem, statistics = "expectedInfluence", order = "sample") + title_theme + ggtitle("Expected Influence – 8w (REM)"),
       width = 8, height = 8, units = "in")

ggsave(filename = file.path(save_path_stability, "EI_0w_NonREM.png"),
       plot = plot(casedrop_boot_pre_nonrem, statistics = "expectedInfluence", order = "sample") + title_theme + ggtitle("Expected Influence – 0w (Non-REM)"),
       width = 8, height = 8, units = "in")

ggsave(filename = file.path(save_path_stability, "EI_8w_NonREM.png"),
       plot = plot(casedrop_boot_post_nonrem, statistics = "expectedInfluence", order = "sample") + title_theme + ggtitle("Expected Influence – 8w (Non-REM)"),
       width = 8, height = 8, units = "in")

# CS-coefficient
get_corStability_text <- function(obj, label) {
  output <- capture.output(corStability(obj, cor = 0.7, statistics = "expectedInfluence", verbose = TRUE))
  return(list(label = label, output = output))
}

results_cs <- list(
  get_corStability_text(casedrop_boot_pre, "Pre All (373)"),
  get_corStability_text(casedrop_boot_post, "Post All (373)"),
  get_corStability_text(casedrop_boot_pre_rem, "Pre REM"),
  get_corStability_text(casedrop_boot_post_rem, "Post REM"),
  get_corStability_text(casedrop_boot_pre_nonrem, "Pre NonREM"),
  get_corStability_text(casedrop_boot_post_nonrem, "Post NonREM")
)

doc <- read_docx()
for (res in results_cs) {
  doc <- doc %>% body_add_par(res$label, style = "heading 2")
  for (line in res$output) {
    doc <- body_add_par(doc, line, style = "Normal")
  }
  doc <- body_add_par(doc, "", style = "Normal")
}
print(doc, target = file.path(save_path_stability, "cs_values.docx"))
cat("CS-coefficient Word file saved.\n")

# Nonparametric bootstrap
set.seed(1234)
nonpar_boot_pre <- bootnet(
  data = df_pre,
  default = "none",
  fun = myWeightedGlasso_bootnet,
  nBoots = 5000,
  gamma = 0.5,
  nCores = 1,
  type = "nonparametric",
  maxErrors = 10,
  errorHandling = "continue",
  weighted = TRUE,
  signed = TRUE,
  verbose = TRUE
)

set.seed(1234)
nonpar_boot_post <- bootnet(
  data = df_post,
  default = "none",
  fun = myWeightedGlasso_bootnet,
  nBoots = 5000,
  gamma = 0.5,
  nCores = 1,
  type = "nonparametric",
  maxErrors = 10,
  errorHandling = "continue",
  weighted = TRUE,
  signed = TRUE,
  verbose = TRUE
)

set.seed(1234)
nonpar_boot_pre_rem <- bootnet(
  data = df_pre_rem,
  default = "none",
  fun = myWeightedGlasso_bootnet_rem,
  nBoots = 5000,
  gamma = 0.5,
  nCores = 1,
  type = "nonparametric",
  maxErrors = 10,
  errorHandling = "continue",
  weighted = TRUE,
  signed = TRUE,
  verbose = TRUE
)

set.seed(1234)
nonpar_boot_pre_nonrem <- bootnet(
  data = df_pre_nonrem,
  default = "none",
  fun = myWeightedGlasso_bootnet_nonrem,
  nBoots = 5000,
  gamma = 0.5,
  nCores = 1,
  type = "nonparametric",
  maxErrors = 10,
  errorHandling = "continue",
  weighted = TRUE,
  signed = TRUE,
  verbose = TRUE
)

set.seed(1234)
nonpar_boot_post_rem <- bootnet(
  data = df_post_rem,
  default = "none",
  fun = myWeightedGlasso_bootnet_rem,
  nBoots = 5000,
  gamma = 0.5,
  nCores = 1,
  type = "nonparametric",
  maxErrors = 10,
  errorHandling = "continue",
  weighted = TRUE,
  signed = TRUE,
  verbose = TRUE
)

set.seed(1234)
nonpar_boot_post_nonrem <- bootnet(
  data = df_post_nonrem,
  default = "none",
  fun = myWeightedGlasso_bootnet_nonrem,
  nBoots = 5000,
  gamma = 0.5,
  nCores = 1,
  type = "nonparametric",
  maxErrors = 10,
  errorHandling = "continue",
  weighted = TRUE,
  signed = TRUE,
  verbose = TRUE
)

# Save confidence interval plots
save_boot_plot <- function(obj, name, title = NULL) {
  custom_title <- ifelse(is.null(title), name, title)
  file_name <- file.path(save_path_accuracy, paste0(name, ".png"))
  p <- plot(obj, order = "sample", plot = "area", verbose = TRUE)
  ggsave(filename = file_name, plot = p, width = 6, height = 8, units = "in", dpi = 300)
  message("Saved:", file_name)
}

save_boot_plot(nonpar_boot_pre, "confidence_intervals_0")
save_boot_plot(nonpar_boot_post, "confidence_intervals_8")
save_boot_plot(nonpar_boot_pre_rem, "confidence_intervals_0_rem")
save_boot_plot(nonpar_boot_pre_nonrem, "confidence_intervals_0_nonrem")
save_boot_plot(nonpar_boot_post_rem, "confidence_intervals_8_rem")
save_boot_plot(nonpar_boot_post_nonrem, "confidence_intervals_8_nonrem")

# Save difference plots
save_diff_plot <- function(obj, what, plot_type, onlyNonZero, name) {
  file_name <- file.path(save_path_diff_test, paste0(name, ".png"))
  p <- plot(obj, what, plot = plot_type, onlyNonZero = onlyNonZero, order = "sample", labels = TRUE, verbose = TRUE)
  ggsave(filename = file_name, plot = p, width = 8, height = 8, units = "in", dpi = 300)
  message("Saved:", file_name)
}

# Edge differences
save_diff_plot(nonpar_boot_pre, "edge", "difference", TRUE, "0_edge_diff")
save_diff_plot(nonpar_boot_post, "edge", "difference", TRUE, "8_edge_diff")
save_diff_plot(nonpar_boot_pre_rem, "edge", "difference", TRUE, "0_rem_edge_diff")
save_diff_plot(nonpar_boot_post_rem, "edge", "difference", TRUE, "8_rem_edge_diff")
save_diff_plot(nonpar_boot_pre_nonrem, "edge", "difference", TRUE, "0_nonrem_edge_diff")
save_diff_plot(nonpar_boot_post_nonrem, "edge", "difference", TRUE, "8_nonrem_edge_diff")

# Strength differences
save_diff_plot(nonpar_boot_pre, "strength", "difference", FALSE, "0_strength_diff")
save_diff_plot(nonpar_boot_post, "strength", "difference", FALSE, "8_strength_diff")
save_diff_plot(nonpar_boot_pre_rem, "strength", "difference", FALSE, "0_rem_strength_diff")
save_diff_plot(nonpar_boot_post_rem, "strength", "difference", FALSE, "8_rem_strength_diff")
save_diff_plot(nonpar_boot_pre_nonrem, "strength", "difference", FALSE, "0_nonrem_strength_diff")
save_diff_plot(nonpar_boot_post_nonrem, "strength", "difference", FALSE, "8_nonrem_strength_diff")

cat("All analyses completed successfully.\n")