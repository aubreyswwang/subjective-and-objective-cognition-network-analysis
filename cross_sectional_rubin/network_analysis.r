## Main network analysis script for paper reproducibility
## Input: 'data' (data frame), 'save_path' (optional character)

# 1. VALIDATE INPUT
if (!exists("data") || is.null(data)) {
  stop("Input dataset 'data' must be defined before running this script.")
}

# 2. SET OUTPUT PATH
if (!exists("save_path")) {
  save_path <- "results"
}
dir.create(save_path, showWarnings = FALSE, recursive = TRUE)

# 3. PREPARE DATA SUBSETS
network_data <- data
network_data <- network_data[network_data$Retained == "retained" & !is.na(network_data$Retained), ]

dataqids_biomarkers_pre  <- network_data[, 13:22]
dataqids_biomarkers_post <- network_data[, 36:45]

# Subgroup data (subgroup == 0: remission; subgroup == 1: non-remission)
subgroup_rem_pre      <- network_data[network_data$subgroup == 0 & !is.na(network_data$subgroup), 13:22]
subgroup_nonrem_pre   <- network_data[network_data$subgroup == 1 & !is.na(network_data$subgroup), 13:22]
subgroup_rem_post     <- network_data[network_data$subgroup == 0 & !is.na(network_data$subgroup), 36:45]
subgroup_nonrem_post  <- network_data[network_data$subgroup == 1 & !is.na(network_data$subgroup), 36:45]

# Extract inverse probability weights (IPW)
weights_global   <- network_data$IPW_w[network_data$Retained == "retained" & !is.na(network_data$Retained)]
weights_rem      <- network_data$IPW_w[network_data$subgroup == 0 & !is.na(network_data$subgroup)]
weights_nonrem   <- network_data$IPW_w[network_data$subgroup == 1 & !is.na(network_data$subgroup)]

# 4. DEFINE VARIABLE NAMES (must be pre-defined in global environment)
# Expected global variables: colnames_abb, groups_vec, groups, node_colors

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

list2env(data_list_pre,  envir = .GlobalEnv)
list2env(data_list_post, envir = .GlobalEnv)

# 5. SELECT AND CLEAN VARIABLES
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

# 6. ESTIMATE MGM NETWORKS
estimate_mgm_network <- function(data) {
  p <- ncol(data)
  fit <- mgm(
    data = data,
    type = rep("g", p),
    level = rep(1, p),
    lambdaSel = "CV",
    ruleReg = "OR"
  )
  pred <- predict(fit, data = data, errorCon = "R2")
  return(list(fit = fit, pred = pred))
}

# Estimate all networks
net_all_pre      <- estimate_mgm_network(biomarkers_qids_pre)
net_all_post     <- estimate_mgm_network(biomarkers_qids_post)
net_rem_pre      <- estimate_mgm_network(biomarkers_rem_pre)
net_nonrem_pre   <- estimate_mgm_network(biomarkers_nonrem_pre)
net_rem_post     <- estimate_mgm_network(biomarkers_rem_post)
net_nonrem_post  <- estimate_mgm_network(biomarkers_nonrem_post)

# 7. CALCULATE MEAN PREDICTABILITY (R²)
mean_predict_biomarkers_qids_pre      <- round(mean(net_all_pre$pred$error$R2), digits = 2)
mean_predict_biomarkers_qids_post     <- round(mean(net_all_post$pred$error$R2), digits = 2)
mean_predict_biomarkers_rem_pre       <- round(mean(net_rem_pre$pred$error$R2), digits = 2)
mean_predict_biomarkers_nonrem_pre    <- round(mean(net_nonrem_pre$pred$error$R2), digits = 2)
mean_predict_biomarkers_rem_post      <- round(mean(net_rem_post$pred$error$R2), digits = 2)
mean_predict_biomarkers_nonrem_post   <- round(mean(net_nonrem_post$pred$error$R2), digits = 2)

# 8. ESTIMATE WEIGHTED EBICGLASSO NETWORKS
compute_weighted_cor <- function(data, weights) {
  cov_mat <- cov.wt(data, wt = weights, cor = FALSE)$cov
  cov2cor(cov_mat)
}

ess_from_weights <- function(w) sum(w)^2 / sum(w^2)

# Compute correlation matrices and effective sample sizes
cormatrix_all_pre    <- compute_weighted_cor(biomarkers_qids_pre, weights_global)
cormatrix_all_post   <- compute_weighted_cor(biomarkers_qids_post, weights_global)
cormatrix_rem_pre    <- compute_weighted_cor(biomarkers_rem_pre, weights_rem)
cormatrix_nonrem_pre <- compute_weighted_cor(biomarkers_nonrem_pre, weights_nonrem)
cormatrix_rem_post   <- compute_weighted_cor(biomarkers_rem_post, weights_rem)
cormatrix_nonrem_post<- compute_weighted_cor(biomarkers_nonrem_post, weights_nonrem)

ess_all     <- ess_from_weights(weights_global)
ess_rem     <- ess_from_weights(weights_rem)
ess_nonrem  <- ess_from_weights(weights_nonrem)

# Network plotting function
plot_qgraph <- function(cor_matrix, predictability, labels, groups_vec, groups, node_colors, sample_size, title = "") {
  qgraph(
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
    vsize = 8,
    edge.labels = TRUE,
    edge.label.cex = 0.65,
    fade = TRUE,
    title = title
  )
}

# Generate all network objects
net_0 <- plot_qgraph(
  cor_matrix = cormatrix_all_pre,
  predictability = net_all_pre$pred,
  labels = colnames_abb,
  groups_vec = groups_vec,
  groups = groups,
  node_colors = node_colors,
  sample_size = ess_all,
  title = "Week 0 (Weighted)"
)

net_8 <- plot_qgraph(
  cor_matrix = cormatrix_all_post,
  predictability = net_all_post$pred,
  labels = colnames_abb,
  groups_vec = groups_vec,
  groups = groups,
  node_colors = node_colors,
  sample_size = ess_all,
  title = "Week 8 (Weighted)"
)

net_rem_0 <- plot_qgraph(
  cor_matrix = cormatrix_rem_pre,
  predictability = net_rem_pre$pred,
  labels = colnames_abb,
  groups_vec = groups_vec,
  groups = groups,
  node_colors = node_colors,
  sample_size = ess_rem,
  title = "Week 0_Remission (Weighted)"
)

net_nonrem_0 <- plot_qgraph(
  cor_matrix = cormatrix_nonrem_pre,
  predictability = net_nonrem_pre$pred,
  labels = colnames_abb,
  groups_vec = groups_vec,
  groups = groups,
  node_colors = node_colors,
  sample_size = ess_nonrem,
  title = "Week 0_Non-remission (Weighted)"
)

net_rem_8 <- plot_qgraph(
  cor_matrix = cormatrix_rem_post,
  predictability = net_rem_post$pred,
  labels = colnames_abb,
  groups_vec = groups_vec,
  groups = groups,
  node_colors = node_colors,
  sample_size = ess_rem,
  title = "Week 8_Remission (Weighted)"
)

net_nonrem_8 <- plot_qgraph(
  cor_matrix = cormatrix_nonrem_post,
  predictability = net_nonrem_post$pred,
  labels = colnames_abb,
  groups_vec = groups_vec,
  groups = groups,
  node_colors = node_colors,
  sample_size = ess_nonrem,
  title = "Week 8_Non-remission (Weighted)"
)

# 9. CENTRALITY ESTIMATES
save_centrality_plot <- function(graph_list, labels, include) {
  centralityPlot(
    graph_list,
    weighted = TRUE,
    signed = TRUE,
    scale = "z-scores",
    labels = labels,
    include = include,
    orderBy = "default"
  )
}

# Expected Influence
centrality_biomarkers_qids_EI <- save_centrality_plot(
  list("Week 0" = net_0, "Week 8" = net_8),
  labels = colnames_abb,
  include = "ExpectedInfluence"
)

centrality_biomarkers_qids_EI_rem <- save_centrality_plot(
  list("Week 0_Remission" = net_rem_0, "Week 8_Remission" = net_rem_8),
  labels = colnames_abb,
  include = "ExpectedInfluence"
)

centrality_biomarkers_qids_EI_nonrem <- save_centrality_plot(
  list("Week 0_Non-remission" = net_nonrem_0, "Week 8_Non-remission" = net_nonrem_8),
  labels = colnames_abb,
  include = "ExpectedInfluence"
)

centrality_biomarkers_qids_EI_week0_comp <- save_centrality_plot(
  list("Week 0_Remission" = net_rem_0, "Week 0_Non-remission" = net_nonrem_0),
  labels = colnames_abb,
  include = "ExpectedInfluence"
)

centrality_biomarkers_qids_EI_week8_comp <- save_centrality_plot(
  list("Week 8_Remission" = net_rem_8, "Week 8_Non-remission" = net_nonrem_8),
  labels = colnames_abb,
  include = "ExpectedInfluence"
)

# Strength
centrality_biomarkers_qids_strength <- save_centrality_plot(
  list("Week 0" = net_0, "Week 8" = net_8),
  labels = colnames_abb,
  include = "Strength"
)

centrality_biomarkers_qids_strength_rem <- save_centrality_plot(
  list("Week 0_Remission" = net_rem_0, "Week 8_Remission" = net_rem_8),
  labels = colnames_abb,
  include = "Strength"
)

centrality_biomarkers_qids_strength_nonrem <- save_centrality_plot(
  list("Week 0_Non-remission" = net_nonrem_0, "Week 8_Non-remission" = net_nonrem_8),
  labels = colnames_abb,
  include = "Strength"
)

centrality_biomarkers_qids_strength_week0_comp <- save_centrality_plot(
  list("Week 0_Remission" = net_rem_0, "Week 0_Non-remission" = net_nonrem_0),
  labels = colnames_abb,
  include = "Strength"
)

centrality_biomarkers_qids_strength_week8_comp <- save_centrality_plot(
  list("Week 8_Remission" = net_rem_8, "Week 8_Non-remission" = net_nonrem_8),
  labels = colnames_abb,
  include = "Strength"
)

# 10. BOOTSTRAP FOR STABILITY (Case-dropping)
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

# Custom bootnet estimators
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

# Run case-dropping bootstrap
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
