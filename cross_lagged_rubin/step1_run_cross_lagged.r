## Batch runner for cross-lagged panel network analysis across multiple imputed datasets
## This script calls cross_lagged_analysis_function.R for each imputed dataset (imp1–imp10).

# 1. SET CRAN MIRROR (for faster package installation in China)
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

# 2. LOAD REQUIRED PACKAGES
pkgs <- c("readxl", "dplyr", "haven", "psych", "car", "NetworkComparisonTest",
          "ggplot2", "mgm", "mice", "bootnet", "qgraph", "psychTools", 
          "glmnet", "lavaan", "readr", "EstimateGroupNetwork", "networktools", 
          "Hmisc", "wCorr", "writexl", "patchwork", "openxlsx", "huge", "extrafont",
          "officer", "flextable", "Matrix", "furrr", "future", "NetworkToolbox")

invisible(lapply(pkgs, function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}))

# 3. SET PATHS (use relative paths for open-source)
data_file      <- "imputed_datasets_all.xlsx"
analysis_script <- "cross_lagged_analysis_function.R"
output_file     <- "cross_lagged_network_results.RData"

# 4. VALIDATE INPUT FILES
if (!file.exists(data_file)) {
  stop(" Input file not found: ", data_file)
}
if (!file.exists(analysis_script)) {
  stop(" Analysis script not found: ", analysis_script)
}

# 5. GET TARGET SHEETS (imp1 to imp10)
sheet_names <- excel_sheets(data_file)
target_sheet_names <- sheet_names[sheet_names %in% paste0("imp", 1:10)]

if (length(target_sheet_names) == 0) {
  stop(" No sheets matching 'imp1' to 'imp10' found in the Excel file.")
}

cat(" Processing sheets:", paste(target_sheet_names, collapse = ", "), "\n")

# 6. RUN ANALYSIS ON EACH SHEET
edges_list <- list()
outEI_list <- list()
inEI_list <- list()
cs_list <- list()
successful_sheets <- c()

for (sheet in target_sheet_names) {
  cat(" Processing sheet:", sheet, "\n")
  
  # Make sheet name available to the analysis function
  assign("sheet_name", sheet, envir = .GlobalEnv)
  
  # Run analysis
  result <- tryCatch({
    source(analysis_script)
  }, error = function(e) {
    cat(" Error in sheet", sheet, ":", conditionMessage(e), "\n")
    return(NULL)
  })
  
  if (!is.null(result)) {
    edges_list[[sheet]] <- result$edges_df
    outEI_list[[sheet]] <- result$outEI_df
    inEI_list[[sheet]] <- result$inEI_df
    cs_list[[sheet]] <- result$cs_df
    successful_sheets <- c(successful_sheets, sheet)
    cat(" Sheet", sheet, "processed successfully.\n")
  } else {
    cat(" Sheet", sheet, "failed.\n")
  }
}

# 7. SAVE RESULTS
if (length(successful_sheets) == 0) {
  stop(" All sheets failed to process.")
}

dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
save(edges_list, outEI_list, inEI_list, cs_list, file = output_file)

cat(" Successfully processed", length(successful_sheets), "sheets:", 
    paste(successful_sheets, collapse = ", "), "\n")
cat(" Results saved to:", output_file, "\n")