# ============================================================
# Assignment 2 — DGE classification
# Files needed (place in ./Raw_Data/):
#   - DEGs_Data_1.csv
#   - DEGs_Data_2.csv
# Each file must have columns: gene_id, logFC, padj
# Output: ./Results/Processed_DEGs_Data_*.csv + console summaries
# ============================================================

# ---------------------------
# 1) Setup: folders & files
# ---------------------------
input_dir  <- "Raw_Data"
output_dir <- "Results"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

files_to_process <- c("DEGs_Data_1.csv", "DEGs_Data_2.csv")

# ---------------------------
# 2) Classification function
# ---------------------------
# Rules:
# - Upregulated   : logFC > 1  & padj < 0.05
# - Downregulated : logFC < -1 & padj < 0.05
# - Otherwise     : Not_Significant
classify_gene <- function(logFC, padj) {
  # Handle missing padj as non-significant
  if (is.na(padj)) padj <- 1
  
  if (logFC > 1 && padj < 0.05) {
    return("Upregulated")
  } else if (logFC < -1 && padj < 0.05) {
    return("Downregulated")
  } else {
    return("Not_Significant")
  }
}

# ---------------------------
# 3) Processing loop
# ---------------------------
result_list <- list()

for (fname in files_to_process) {
  cat("\n============================\n")
  cat("Processing:", fname, "\n")
  
  # Build full path and read
  fpath <- file.path(input_dir, fname)
  if (!file.exists(fpath)) {
    stop(paste("File not found:", fpath,
               "\nMake sure it is inside the Raw_Data folder."))
  }
  
  # Read CSV (no factors to avoid surprises)
  dat <- read.csv(fpath, header = TRUE, stringsAsFactors = FALSE)
  
  # --- Basic checks on required columns ---
  required_cols <- c("gene_id", "logFC", "padj")
  missing_cols  <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop(paste("The file", fname, "is missing required column(s):",
               paste(missing_cols, collapse = ", ")))
  }
  
  # --- Handle missing padj as 1 (not significant) ---
  na_padj_before <- sum(is.na(dat$padj))
  if (na_padj_before > 0) {
    dat$padj[is.na(dat$padj)] <- 1
    cat("Replaced", na_padj_before, "missing padj values with 1.\n")
  }
  
  # --- Add status column using classify_gene ---
  dat$status <- mapply(classify_gene, dat$logFC, dat$padj)
  
  # --- Print summaries ---
  cat("Status counts (Up/Down/Not_Significant):\n")
  print(table(dat$status))
  
  # Significant = padj < 0.05 (regardless of direction)
  dat$significant <- dat$padj < 0.05
  cat("Significant (padj < 0.05) count:\n")
  print(table(dat$significant))
  
  # Optional: Up/Down among significant only
  up_count   <- sum(dat$status == "Upregulated")
  down_count <- sum(dat$status == "Downregulated")
  sig_count  <- sum(dat$significant)
  cat("Summary → Significant:", sig_count,
      "| Upregulated:", up_count,
      "| Downregulated:", down_count, "\n")
  
  # --- Save processed file ---
  outpath <- file.path(output_dir, paste0("Processed_", fname))
  write.csv(dat, outpath, row.names = FALSE)
  cat("Saved:", outpath, "\n")
  
  # Keep in memory too
  result_list[[fname]] <- dat
}

# ---------------------------
# 4) Access in R if needed
# ---------------------------
# results_1 <- result_list[["DEGs_Data_1.csv"]]
# results_2 <- result_list[["DEGs_Data_2.csv"]]

cat("\nAll done. Check the Results folder for processed files.\n")
