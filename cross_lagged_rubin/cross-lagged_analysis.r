## Cross-lagged panel network analysis for a single imputed dataset
## This function is designed to be called within a loop over multiple imputations.

# Input: sheet_name (character) - name of the sheet in imputed_datasets_all.xlsx
# Output: list containing edge weights, centrality estimates, and stability coefficients

# 1. LOAD DATA
CoV_l_2 <- read_excel("imputed_datasets_all.xlsx", sheet = sheet_name)
CoV_l_2 <- CoV_l_2 %>%
  filter(Retained == "retained" & !is.na(Retained))

weights_global <- CoV_l_2$IPW_w

# 2. PREPARE DATA
w0b <- as.data.frame(CoV_l_2[, 13:22])
w8b <- as.data.frame(CoV_l_2[, 36:45])

colnames(w0b) <- colnames(w8b) <- c(
  "S-AT", "S-ME", "S-EF",
  "O-ME", "O-EF", "O-AT", "O-PS",
  "SM", "GI", "EN"
)

w0b <- as.data.frame(lapply(w0b, as.numeric))
w8b <- as.data.frame(lapply(w8b, as.numeric))

# Combine with covariates
w0_8 <- data.frame(w0b, w8b)
w0_8$age <- CoV_l_2$Age
w0_8 <- huge.npn(w0_8)
w0_8$gender <- as.numeric(CoV_l_2$Sex)

# 3. ESTIMATE CROSS-LAGGED NETWORK USING GLMNET
numCovar <- 2  # age + gender
k <- 10         # number of nodes

adjMat3 <- matrix(0, nrow = k + numCovar, ncol = k + numCovar)
w0_8 <- makeX(as.data.frame(w0_8))
x_matrix <- w0_8[, c(1:k, (2*k + 1):(2*k + numCovar))]

for (i in 1:k) {
  set.seed(1)
  fit <- cv.glmnet(
    x           = x_matrix,
    y           = w0_8[, k + i],
    family      = "gaussian",
    alpha       = 1,
    standardize = TRUE,
    weights     = weights_global
  )
  lam <- fit$lambda.min
  adjMat3[, i] <- coef(fit, s = lam)[-1]  # exclude intercept
}

# Extract T1 -> T2 adjacency matrix
labels <- c(colnames(w0b), "age", "gender")
adjMat3_full <- getWmat(adjMat3, nNodes = k + numCovar, labels = labels, directed = TRUE)
adjMat3_T1T2 <- adjMat3_full[1:k, 1:k]

# Remove self-loops
diag(adjMat3_T1T2) <- 0

# 4. CENTRALITY ESTIMATION
groups_vec <- c(
  rep("Subjective Cognition", 3),
  rep("Objective Cognition", 4),
  rep("Depressive Symptom", 3)
)

nwc2 <- qgraph(
  adjMat3_T1T2,
  groups = groups_vec,
  legend = FALSE,
  threshold = 0.05,
  labels = colnames(w0b)
)

plot_in <- centralityPlot(
  list("Week 0" = nwc2),
  weighted = TRUE,
  signed = TRUE,
  scale = "z-scores",
  labels = colnames(w0b),
  include = "InExpectedInfluence",
  orderBy = "InExpectedInfluence"
)

plot_out <- centralityPlot(
  list("Week 0" = nwc2),
  weighted = TRUE,
  signed = TRUE,
  scale = "z-scores",
  labels = colnames(w0b),
  include = "OutExpectedInfluence",
  orderBy = "OutExpectedInfluence"
)

# 5. BOOTSTRAP FOR STABILITY
CLPN.OR <- function(data) {
  k <- ncol(w0b)
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
    adj[, i] <- coef(fit, s = lam)[-1]  # exclude intercept
  }
  return(adj)
}

boot_data <- data.frame(w0b, w8b, ipw = weights_global)
boot_data <- makeX(boot_data)
net.2 <- estimateNetwork(boot_data, fun = CLPN.OR, labels = colnames(w0b), directed = TRUE)

caseBoot.w2 <- bootnet(
  net.2,
  type = "case",
  nBoots = 1000,
  directed = TRUE,
  statistics = c("outExpectedInfluence", "inExpectedInfluence")
)

# 6. EXTRACT RESULTS
outEI_values <- plot_out$data
inEI_values <- plot_in$data

cs_in  <- corStability(caseBoot.w2, statistics = "inExpectedInfluence")
cs_out <- corStability(caseBoot.w2, statistics = "outExpectedInfluence")

# Edge list
diag(adjMat3_T1T2) <- 1  # temporarily restore for ranking
res.w2 <- order(adjMat3_T1T2, decreasing = TRUE)[seq_len(400)]
pos.w2 <- arrayInd(res.w2, dim(adjMat3_T1T2), useNames = TRUE)
posWithLabs.w2 <- data.frame(
  nodeOut   = colnames(w0b)[pos.w2[, 1]],
  nodeIn    = colnames(w0b)[pos.w2[, 2]],
  value     = adjMat3_T1T2[res.w2],
  DomainOut = groups_vec[pos.w2[, 1]],
  DomainIn  = groups_vec[pos.w2[, 2]]
)
all_edges <- data.frame(posWithLabs.w2)

# 7. CONSTRUCT OUTPUT
edges_df <- data.frame(
  edge = all_edges,
  EdgeWeight = as.vector(adjMat3_T1T2)
)

outEI_df <- data.frame(
  Node = as.character(outEI_values$node),
  OutEI = outEI_values$value
)

inEI_df <- data.frame(
  Node = as.character(inEI_values$node),
  InEI = inEI_values$value
)

cs_df <- data.frame(
  Centrality = c("InExpectedInfluence", "OutExpectedInfluence"),
  Stability_Coefficient = c(cs_in, cs_out)
)

# Return structured results for Rubin merging
list(
  edges_df = edges_df,
  outEI_df = outEI_df,
  inEI_df = inEI_df,
  cs_df = cs_df
)