## Imputation, Inverse Probability Weighting (IPW), and Mean Imputed Dataset

# 1. SET OUTPUT DIRECTORY (relative path)
out_dir <- "imputation"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 2. LOAD REQUIRED PACKAGES
pkgs <- c(
  "readxl", "dplyr", "mice", "VIM", "naniar", "ggplot2", "lattice",
  "grid", "ipw", "cobalt", "MatchThem", "survey", "openxlsx", "qgraph",
  "bootnet", "WeightIt"
)

# Install missing packages
invisible(lapply(pkgs, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}))

# Load all packages
lapply(pkgs, library, character.only = TRUE)

# 3. READ AND PREPROCESS DATA
df <- read_excel("data.xlsx", sheet = "sheet1") %>%
  mutate(
    ID = as.character(ID),
    Cohort = factor(Cohort, levels = c("A", "B")),
    Retained = factor(Retained, levels = c("retained", "lost")),
    Sex = factor(Sex),
    Marriage_Status = factor(Marriage_Status),
    Family_Income = factor(Family_Income),
    Rural = factor(Rural),
    First_Episode = factor(First_Episode),
    Education = as.numeric(Education),
    Disease_Duration = as.numeric(Disease_Duration)
  )

# 4. DEFINE VARIABLE TYPES FOR IMPUTATION
# Ordered categorical variables (baseline)
ord9_w0 <- c("w0_Memory", "w0_Planning_Organisation")
ord5_w0 <- c("w0_Attention_Concentration")
ord4_w0 <- c(
  "w0_Sleep_Disorder", "w0_Sad_Mood", "w0_Appetite_Weight_Change",
  "w0_Concentration_Decision_Making", "w0_Outlook", "w0_Suicidal_Ideation",
  "w0_General_Interest", "w0_Energy_Fatigability", "w0_Psychomotor_Agitation_Retardation"
)

# Ordered categorical variables (week 8)
ord9_w8 <- sub("^w0_", "w8_", ord9_w0)
ord5_w8 <- sub("^w0_", "w8_", ord5_w0)
ord4_w8 <- sub("^w0_", "w8_", ord4_w0)

# Apply ordered factors
df <- df %>%
  mutate(
    across(all_of(ord9_w0), ~ ordered(.x, levels = seq(0, 4, by = 0.5))),
    across(all_of(ord5_w0), ~ ordered(.x, levels = 0:4)),
    across(all_of(ord4_w0), ~ ordered(.x, levels = 0:3)),
    across(all_of(ord9_w8), ~ ordered(.x, levels = seq(0, 4, by = 0.5))),
    across(all_of(ord5_w8), ~ ordered(.x, levels = 0:4)),
    across(all_of(ord4_w8), ~ ordered(.x, levels = 0:3))
  )

# Ensure remaining w0_/w8_ variables are numeric
w0_vars <- grep("^w0_", names(df), value = TRUE)
w8_vars <- grep("^w8_", names(df), value = TRUE)

scale_vars_w0 <- setdiff(w0_vars, c(ord9_w0, ord5_w0, ord4_w0, "w0_HAMD", "w0_MADRS"))
scale_vars_w8 <- setdiff(w8_vars, c(ord9_w8, ord5_w8, ord4_w8, "w8_HAMD", "w8_MADRS"))

df[scale_vars_w0] <- lapply(df[scale_vars_w0], as.numeric)
df[scale_vars_w8] <- lapply(df[scale_vars_w8], as.numeric)

# 5. MISSINGNESS PLOTS
# Helper: Mask cohort-specific missingness in HAMD/MADRS
mask_cohort_missing <- function(data) {
  data %>%
    mutate(
      w0_HAMD  = if_else(Cohort == "B", NA_real_, w0_HAMD),
      w0_MADRS = if_else(Cohort == "A", NA_real_, w0_MADRS),
      w8_HAMD  = if_else(Cohort == "B", NA_real_, w8_HAMD),
      w8_MADRS = if_else(Cohort == "A", NA_real_, w8_MADRS)
    )
}

# Full missingness pattern
png(file.path(out_dir, "missing_pattern_full.png"), width = 800, height = 600, res = 100)
viz_full <- mask_cohort_missing(df)
gg_miss_upset(viz_full %>% select(-Cohort), nsets = 8)
mask_full <- is.na(viz_full)
colnames(mask_full) <- names(viz_full)
mask_full[viz_full$Cohort == "B", c("w0_HAMD", "w8_HAMD")] <- FALSE
mask_full[viz_full$Cohort == "A", c("w0_MADRS", "w8_MADRS")] <- FALSE
miss_pct_full <- colMeans(mask_full[, !colnames(mask_full) %in% c("Cohort", "w0_HAMD", "w0_MADRS", "w8_HAMD", "w8_MADRS")])
par(mar = c(15, 4, 4, 2) + 0.1)
barplot(sort(miss_pct_full, TRUE), las = 2, cex.names = 0.7,
        main = "Within-cohort Missing Rates", ylab = "Missing Proportion")
dev.off()

# Baseline only
png(file.path(out_dir, "missing_pattern_w0.png"), width = 800, height = 600, res = 100)
viz_w0 <- df %>% select(all_of(w0_vars), Cohort) %>% mask_cohort_missing()
gg_miss_upset(viz_w0 %>% select(-Cohort), nsets = 8)
mask_w0 <- is.na(viz_w0)
mask_w0[viz_w0$Cohort == "B", "w0_HAMD"] <- FALSE
mask_w0[viz_w0$Cohort == "A", "w0_MADRS"] <- FALSE
miss_pct_w0 <- colMeans(mask_w0[, !colnames(mask_w0) %in% c("Cohort", "w0_HAMD", "w0_MADRS")])
par(mar = c(15, 4, 4, 2) + 0.1)
barplot(sort(miss_pct_w0, TRUE), las = 2, cex.names = 0.7,
        main = "Baseline Missing Rates", ylab = "Missing Proportion")
dev.off()

# Week 8 only
png(file.path(out_dir, "missing_pattern_w8.png"), width = 800, height = 600, res = 100)
viz_w8 <- df %>% select(all_of(w8_vars), Cohort) %>% mask_cohort_missing()
gg_miss_upset(viz_w8 %>% select(-Cohort), nsets = 8)
mask_w8 <- is.na(viz_w8)
mask_w8[viz_w8$Cohort == "B", "w8_HAMD"] <- FALSE
mask_w8[viz_w8$Cohort == "A", "w8_MADRS"] <- FALSE
miss_pct_w8 <- colMeans(mask_w8[, !colnames(mask_w8) %in% c("Cohort", "w8_HAMD", "w8_MADRS")])
par(mar = c(15, 4, 4, 2) + 0.1)
barplot(sort(miss_pct_w8, TRUE), las = 2, cex.names = 0.7,
        main = "Week 8 Missing Rates", ylab = "Missing Proportion")
dev.off()

# 6. MULTIPLE IMPUTATION SETUP
no_imp <- c(
  "ID", "subgroup", "Retained",
  "w0_HAMD", "w0_MADRS", "w0_Prospective_Memory", "w0_Retrospective_Memory",
  "w0_Inhibitory_Control", "w0_Working_Memory", "w0_Cognitive_Flexibility",
  "w8_HAMD", "w8_MADRS", "w8_Prospective_Memory", "w8_Retrospective_Memory",
  "w8_Inhibitory_Control", "w8_Working_Memory", "w8_Cognitive_Flexibility"
)

ord_vars <- c(ord5_w0, ord4_w0, ord9_w0, ord5_w8, ord4_w8, ord9_w8)
num_vars <- setdiff(names(df)[sapply(df, is.numeric)], c(no_imp, ord_vars))
fact_vars <- c("Sex", "Marriage_Status", "Family_Income", "Rural", "First_Episode")

# Build imputation settings
meth <- make.method(df)
pred <- make.predictorMatrix(df)

meth[no_imp] <- ""
pred[, no_imp] <- 0

# Cohort: predictor only, not imputed
meth["Cohort"] <- ""
pred["Cohort", ] <- 0
pred[, "Cohort"] <- 1

# Imputation methods
meth[num_vars] <- "pmm"
meth[ord_vars] <- "polr"
meth[fact_vars] <- "logreg"
multilevel <- fact_vars[sapply(df[fact_vars], nlevels) > 2]
meth[multilevel] <- "polyreg"

# Run imputation
set.seed(123)
imp <- mice(
  df,
  method = meth,
  predictorMatrix = pred,
  where = is.na(df) & !names(df) %in% no_imp,
  m = 10,
  maxit = 10,
  printFlag = TRUE
)

# Diagnostic plots
plot(imp, c("Age", "Disease_Duration"))
for (v in num_vars) {
  png(file.path(out_dir, paste0("density_", v, ".png")), width = 1200, height = 800, res = 100)
  print(densityplot(imp, as.formula(paste0("~", v)),
                    main = list(paste0(v, " density"), cex = 2.5),
                    xlab = list(v, cex = 2),
                    ylab = list("Density", cex = 2),
                    scales = list(cex = 1.5),
                    plot.points = FALSE))
  dev.off()
}

# 7. INVERSE PROBABILITY WEIGHTING (IPW) FOR ATTRITION BIAS
ps_weights <- weightthem(
  Retained ~ Age + Sex + Marriage_Status + Rural +
             Education + Disease_Duration + First_Episode + Family_Income,
  by = ~ Cohort,
  data = imp,
  method = "ps",
  approach = "within",
  estimand = "ATE"
)

bal_tab <- bal.tab(ps_weights, un = TRUE)
balance_df <- bal_tab$Balance.Across.Imputations
balance_df <- balance_df[!grepl("prop.score", rownames(balance_df), ignore.case = TRUE), ]

covariates <- rownames(balance_df)
unadjusted <- balance_df$Mean.Diff.Un
adjusted <- balance_df$Mean.Diff.Adj

# Love plot
png(file.path(out_dir, "balance_loveplot.png"), width = 1200, height = 800, res = 100)
par(mar = c(5, 12, 4, 2))
max_diff <- max(abs(c(unadjusted, adjusted)), na.rm = TRUE)
xlim_range <- c(-max_diff - 0.1, max_diff + 0.1)
plot(1, 1, type = "n", xlim = xlim_range, ylim = c(0.5, length(covariates) + 0.5),
     xlab = "Standardized Mean Difference", ylab = "",
     main = "Covariate Balance Before and After Weighting",
     yaxt = "n", cex.main = 2.5, cex.lab = 2, cex.axis = 1.5)
abline(v = c(-0.1, 0, 0.1), lty = c(2, 1, 2), col = c("red", "black", "red"))
abline(v = seq(-1, 1, 0.05), lty = 3, col = "lightgray")
points(unadjusted, 1:length(covariates), pch = 1, col = "red", cex = 2)
points(adjusted, 1:length(covariates), pch = 16, col = "blue", cex = 2)
axis(2, at = 1:length(covariates), labels = covariates, las = 1, cex.axis = 1.5)
legend("topright", legend = c("Unweighted", "Weighted"),
       col = c("red", "blue"), pch = c(1, 16), bty = "n", cex = 1.5)
dev.off()

# Export SMD table
bal_mat <- as.data.frame(bal_tab$Balance.Across.Imputations)
bal_mat <- bal_mat[!grepl("prop.score", rownames(bal_mat), ignore.case = TRUE), ]
bal_mat$Variable <- rownames(bal_mat)
diff_cols <- grep("Diff", names(bal_mat), value = TRUE)
bal_df <- bal_mat %>% select(
  Variable,
  SMD_Unadjusted = all_of(diff_cols[2]),
  SMD_Adjusted   = all_of(diff_cols[5])
)
write.csv(bal_df, file.path(out_dir, "balance_smd.csv"), row.names = FALSE)

# 8. EXTRACT AND STABILIZE IPW WEIGHTS
overall_retained_prob <- mean(df$Retained == "retained", na.rm = TRUE)

wt_list <- lapply(seq_len(imp$m), function(k) {
  w_obj <- ps_weights$models[[k]]
  df_k  <- complete(imp, k)
  w_raw <- w_obj$weights
  Rnum  <- as.numeric(df_k$Retained == "retained")
  w_stable <- ifelse(Rnum == 1, w_raw * overall_retained_prob, w_raw * (1 - overall_retained_prob))
  q_low <- quantile(w_stable, 0.01, na.rm = TRUE)
  q_up  <- quantile(w_stable, 0.99, na.rm = TRUE)
  w_trunc <- pmax(q_low, pmin(q_up, w_stable))
  data.frame(ID = df_k$ID, IPW_w = w_trunc, stringsAsFactors = FALSE)
})

# Compute average weights
w_long <- bind_rows(wt_list)
avg_wt <- w_long %>%
  group_by(ID) %>%
  summarise(IPW_w = mean(IPW_w, na.rm = TRUE), .groups = "drop")

# 9. CREATE MEAN IMPUTED DATASET
all_imp_df <- lapply(seq_len(imp$m), function(k) complete(imp, k))
mean_imputed_data <- all_imp_df[[1]] %>% select(ID) %>% distinct()

for (var in names(all_imp_df[[1]])) {
  if (var == "ID") next
  if (is.numeric(all_imp_df[[1]][[var]])) {
    mean_val <- sapply(seq_len(nrow(mean_imputed_data)), function(i) {
      mean(sapply(all_imp_df, function(df) df[[var]][df$ID == mean_imputed_data$ID[i]]), na.rm = TRUE)
    })
    mean_imputed_data[[var]] <- mean_val
  } else {
    mode_val <- sapply(seq_len(nrow(mean_imputed_data)), function(i) {
      vals <- unlist(lapply(all_imp_df, function(df) df[[var]][df$ID == mean_imputed_data$ID[i]]))
      names(sort(table(vals), decreasing = TRUE))[1]
    })
    mean_imputed_data[[var]] <- as.factor(mode_val)
  }
}

# Merge average IPW weights
mean_imputed_data <- left_join(mean_imputed_data, avg_wt, by = "ID")

# 10. SAVE RESULTS
# Save mean imputed dataset (for network analysis)
write.xlsx(
  list(Mean_Imputed = mean_imputed_data),
  file = file.path("data", "mean_imputed_with_single_ipw.xlsx")
)

# Save all imputed datasets (optional, for transparency)
imp_datasets <- lapply(seq_len(imp$m), function(k) {
  df_k <- complete(imp, k)
  wt_k <- wt_list[[k]]
  left_join(df_k, wt_k, by = "ID")
})
names(imp_datasets) <- paste0("imp", seq_len(imp$m))
write.xlsx(imp_datasets, file = file.path(out_dir, "imputed_datasets_all.xlsx"))

# Save IPW weight summary
weights_all <- unlist(lapply(wt_list, function(x) x$IPW_w))
weight_summary <- data.frame(
  Mean = mean(weights_all, na.rm = TRUE),
  SD = sd(weights_all, na.rm = TRUE),
  Median = median(weights_all, na.rm = TRUE),
  IQR = IQR(weights_all),
  Min = min(weights_all, na.rm = TRUE),
  Max = max(weights_all, na.rm = TRUE),
  `1%` = quantile(weights_all, 0.01, na.rm = TRUE),
  `99%` = quantile(weights_all, 0.99, na.rm = TRUE)
)
write.csv(weight_summary, file.path(out_dir, "ipw_weight_summary.csv"), row.names = FALSE)

cat("Imputation and IPW completed. Mean imputed dataset saved to: mean_imputed_with_single_ipw.xlsx\n")