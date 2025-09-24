## Main entry script for network analysis across multiple imputed datasets
## This script loads imputed datasets and runs the core network analysis on each.

# 1. SET CRAN MIRROR (for faster package installation in China)
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

# 2. LOAD REQUIRED PACKAGES
pkgs <- c("readxl", "dplyr", "haven", "psych", "car", "NetworkComparisonTest",
          "ggplot2", "mgm", "mice", "bootnet", "qgraph", "psychTools", 
          "glmnet", "lavaan", "readr", "EstimateGroupNetwork", "networktools", 
          "Hmisc", "wCorr", "writexl", "patchwork", "openxlsx", "huge", "extrafont",
          "officer", "flextable", "Matrix", "furrr", "future")

invisible(lapply(pkgs, function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}))
library(NetworkToolbox)

# 3. SET DATA PATHS (use relative paths for open-source)
# Note: Assumes this script is run from project root
folder_path <- "data"  # Base data directory
excel_file  <- "imputed_datasets_all.xlsx"
full_path   <- file.path(folder_path, excel_file)

# Create output directory
output_base <- "network_analysis"
dir.create(output_base, showWarnings = FALSE, recursive = TRUE)

network_script_path <- "network_analysis.R"  # Core analysis script

# 4. LOAD IMPUTED DATASETS
sheet_names <- excel_sheets(full_path)[1:10]  # Use first 10 sheets

list_of_datasets <- lapply(sheet_names, function(sheet) {
  read_excel(full_path, sheet = sheet)
})

# 5. DEFINE GLOBAL VARIABLES FOR NETWORK ANALYSIS
colnames_abb <- c(
  "S-AT", "S-ME", "S-EF",
  "O-ME", "O-EF", "O-AT", "O-PS",
  "SM", "GI", "EN"
)

groups_vec <- c(
  rep("Subjective Cognition", 3),
  rep("Objective Cognition", 4),
  rep("Depressive Symptom", 3)
)

group_colors <- c(
  "Subjective Cognition" = "#8BD000FF",
  "Objective Cognition"  = "#E8C32EFF",
  "Depressive Symptom"   = "#F0C9C0FF"
)

items_clean <- list(
  "Attention/Concentration", "Memory", "Planning/Organisation",
  "Memory", "Executive Functions", "Attention", "Processing Speed",
  "Sad Mood", "General Interest", "Energy/Fatigability"
)

node_colors <- group_colors[groups_vec]
groups      <- split(items_clean, groups_vec)

# 6. ANALYSIS FUNCTION FOR EACH IMPUTED DATASET
analyze_dataset <- function(data, index) {
  
  cat("Processing sheet:", sheet_names[index], "\n")
  
  # Create output directory for this imputation
  output_dir <- file.path(output_base, paste0("sheet_", index))
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Assign global variables for network_analysis.R
  assign("data", data, envir = .GlobalEnv)
  assign("save_path", output_dir, envir = .GlobalEnv)
  assign("colnames_abb", colnames_abb, envir = .GlobalEnv)
  assign("groups_vec", groups_vec, envir = .GlobalEnv)
  assign("group_colors", group_colors, envir = .GlobalEnv)
  assign("node_colors", node_colors, envir = .GlobalEnv)
  
  # Run core network analysis
  tryCatch({
    
    source(network_script_path)
    
    # Extract adjacency matrices
    net_0_adj  <- getWmat(net_0)
    net_8_adj  <- getWmat(net_8)
    
    # Extract Expected Influence
    ei_data_pre  <- centrality_biomarkers_qids_EI$data %>%
      filter(type == "Week 0", measure == "ExpectedInfluence") %>%
      dplyr::select(node, value)
    
    ei_data_post <- centrality_biomarkers_qids_EI$data %>%
      filter(type == "Week 8", measure == "ExpectedInfluence") %>%
      dplyr::select(node, value)
    
    # Extract Strength
    strength_data_pre  <- centrality_biomarkers_qids_strength$data %>%
      filter(type == "Week 0", measure == "Strength") %>%
      dplyr::select(node, value)
    
    strength_data_post <- centrality_biomarkers_qids_strength$data %>%
      filter(type == "Week 8", measure == "Strength") %>%
      dplyr::select(node, value)
    
    # Helper: Extract CS-coefficient
    get_cs_value <- function(bootnet_obj) {
      if (is.null(bootnet_obj)) return(NA)
      cs_result <- tryCatch({
        corStability(bootnet_obj, cor = 0.7, statistics = "expectedInfluence", verbose = FALSE)
      }, error = function(e) NA)
      
      if (is.list(cs_result) && !is.null(cs_result$corStability)) {
        cs_result$corStability[[1]]$Estimates[["CS-coefficient"]]
      } else if (is.numeric(cs_result)) {
        cs_result
      } else if (is.character(cs_result)) {
        match <- grep("CS-coefficient =", cs_result, value = TRUE)
        if (length(match) > 0) {
          as.numeric(sub(".*= ", "", match))
        } else {
          NA
        }
      } else {
        NA
      }
    }
    
    # Extract CS-coefficients
    cs_pre           <- get_cs_value(casedrop_boot_pre)
    cs_post          <- get_cs_value(casedrop_boot_post)
    cs_pre_rem       <- get_cs_value(casedrop_boot_pre_rem)
    cs_post_rem      <- get_cs_value(casedrop_boot_post_rem)
    cs_pre_nonrem    <- get_cs_value(casedrop_boot_pre_nonrem)
    cs_post_nonrem   <- get_cs_value(casedrop_boot_post_nonrem)
    
    # Build result list
    result <- list(
      baseline_all = list(
        edge_matrix     = net_0_adj,
        centrality      = ei_data_pre,
        strength        = strength_data_pre,
        predictability  = predictabil_biomarkers_qids_pre$errors,
        weights         = weights_global,
        ess             = ess_global,
        cs              = cs_pre
      ),
      week8_all = list(
        edge_matrix     = net_8_adj,
        centrality      = ei_data_post,
        strength        = strength_data_post,
        predictability  = predictabil_biomarkers_qids_post$errors,
        weights         = weights_global,
        ess             = ess_global,
        cs              = cs_post
      )
    )
    
    # Add subgroup results if available
    other_nets <- c("net_rem_0", "net_rem_8", "net_nonrem_0", "net_nonrem_8")
    other_preds <- c("predictabil_biomarkers_rem_pre", "predictabil_biomarkers_rem_post",
                     "predictabil_biomarkers_nonrem_pre", "predictabil_biomarkers_nonrem_post")
    
    if (all(other_nets %in% ls(.GlobalEnv)) && all(other_preds %in% ls(.GlobalEnv))) {
      
      net_rem_0_adj      <- getWmat(get("net_rem_0"))
      net_rem_8_adj      <- getWmat(get("net_rem_8"))
      net_nonrem_0_adj   <- getWmat(get("net_nonrem_0"))
      net_nonrem_8_adj   <- getWmat(get("net_nonrem_8"))
      
      ei_rem_pre         <- centrality_biomarkers_qids_EI_rem$data %>% filter(type == "Week 0_Remission", measure == "ExpectedInfluence") %>% dplyr::select(node, value)
      ei_rem_post        <- centrality_biomarkers_qids_EI_rem$data %>% filter(type == "Week 8_Remission", measure == "ExpectedInfluence") %>% dplyr::select(node, value)
      ei_nonrem_pre      <- centrality_biomarkers_qids_EI_nonrem$data %>% filter(type == "Week 0_Non-remission", measure == "ExpectedInfluence") %>% dplyr::select(node, value)
      ei_nonrem_post     <- centrality_biomarkers_qids_EI_nonrem$data %>% filter(type == "Week 8_Non-remission", measure == "ExpectedInfluence") %>% dplyr::select(node, value)
      
      strength_rem_pre   <- centrality_biomarkers_qids_strength_rem$data %>% filter(type == "Week 0_Remission", measure == "Strength") %>% dplyr::select(node, value)
      strength_rem_post  <- centrality_biomarkers_qids_strength_rem$data %>% filter(type == "Week 8_Remission", measure == "Strength") %>% dplyr::select(node, value)
      strength_nonrem_pre<- centrality_biomarkers_qids_strength_nonrem$data %>% filter(type == "Week 0_Non-remission", measure == "Strength") %>% dplyr::select(node, value)
      strength_nonrem_post<-centrality_biomarkers_qids_strength_nonrem$data %>% filter(type == "Week 8_Non-remission", measure == "Strength") %>% dplyr::select(node, value)
      
      result$baseline_rem <- list(
        edge_matrix     = net_rem_0_adj,
        centrality      = ei_rem_pre,
        strength        = strength_rem_pre,
        predictability  = predictabil_biomarkers_rem_pre$errors,
        weights         = weights_rem,
        ess             = ess_rem,
        cs              = cs_pre_rem
      )
      result$week8_rem <- list(
        edge_matrix     = net_rem_8_adj,
        centrality      = ei_rem_post,
        strength        = strength_rem_post,
        predictability  = predictabil_biomarkers_rem_post$errors,
        weights         = weights_rem,
        ess             = ess_rem,
        cs              = cs_post_rem
      )
      result$baseline_nonrem <- list(
        edge_matrix     = net_nonrem_0_adj,
        centrality      = ei_nonrem_pre,
        strength        = strength_nonrem_pre,
        predictability  = predictabil_biomarkers_nonrem_pre$errors,
        weights         = weights_nonrem,
        ess             = ess_nonrem,
        cs              = cs_pre_nonrem
      )
      result$week8_nonrem <- list(
        edge_matrix     = net_nonrem_8_adj,
        centrality      = ei_nonrem_post,
        strength        = strength_nonrem_post,
        predictability  = predictabil_biomarkers_nonrem_post$errors,
        weights         = weights_nonrem,
        ess             = ess_nonrem,
        cs              = cs_post_nonrem
      )
      
    } else {
      cat("Subgroup network variables missing; analyzing full sample only.\n")
    }
    
    # Save intermediate results
    save(result, file = file.path(output_dir, "results.RData"))
    cat("Sheet", sheet_names[index], "processed successfully.\n\n")
    return(result)
    
  }, error = function(e) {
    cat("Sheet", sheet_names[index], "failed:", conditionMessage(e), "\n")
    return(NULL)
  })
}

# 7. RUN ANALYSIS ON ALL IMPUTED DATASETS
results_per_dataset <- list()
successful_sheets <- c()

for (i in seq_along(list_of_datasets)) {
  result <- analyze_dataset(list_of_datasets[[i]], i)
  if (!is.null(result)) {
    results_per_dataset[[length(results_per_dataset) + 1]] <- result
    successful_sheets <- c(successful_sheets, sheet_names[i])
  }
}

if (length(results_per_dataset) == 0) {
  stop("All datasets failed to process.")
}

names(results_per_dataset) <- successful_sheets

# 8. SAVE FINAL RESULTS
save(results_per_dataset, file = file.path(output_base, "results_per_dataset.RData"))
cat("All results saved to:", file.path(output_base, "results_per_dataset.RData"), "\n")