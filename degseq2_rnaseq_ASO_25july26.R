# ==========================================
# STEP 1: Dependencies Check & Loading
# ==========================================
cat("Checking and installing packages if missing...\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery")
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")

library(GEOquery)
library(tidyverse)
library(DESeq2)

gse_id <- "GSE270454" 

# ==========================================
# STEP 2: Fetch and Map Metadata
# ==========================================
cat("Fetching and cleaning metadata...\n")

options(timeout = 600)
options(download.file.method.GEOquery = "libcurl")

metadata <- tryCatch({
  gse <- getGEO(gse_id, GSEMatrix = TRUE, getGPL = FALSE)
  pData(gse[])
}, error = function(e) {
  cat("Online connection failed. Attempting offline/memory file lookup...\n")
  local_path <- "GSE270454/GSE270454_series_matrix.txt.gz"
  local_root <- "GSE270454_series_matrix.txt.gz"
  
  if (file.exists(local_path)) {
    pData(getGEO(filename = local_path, getGPL = FALSE))
  } else if (file.exists(local_root)) {
    pData(getGEO(filename = local_root, getGPL = FALSE))
  } else if (exists("metadata")) {
    cat("Success: Restored metadata matrix straight from active R memory.\n")
    return(metadata)
  } else {
    stop("Could not find metadata! Please verify your active session memory or ensure 'GSE270454_series_matrix.txt.gz' is in your directory.")
  }
})

metadata_rna <- metadata %>%
  filter(library_strategy == "RNA-Seq") %>%
  mutate(
    sample_num = str_extract(title, "[0-9]+$"),
    short_prefix = case_when(
      str_detect(title, "middle-aged/older adult") ~ "ASM",
      str_detect(title, "older adult")             ~ "ASO",
      str_detect(title, "MCI patient")             ~ "MCI",
      str_detect(title, "AD patient")              ~ "AD",
      TRUE                                         ~ "Unknown"
    ),
    count_col_id = paste0(short_prefix, "_", sample_num),
    # CHANGES MADE HERE: Only assign ASO as CTL. Explicitly tag ASM to drop it.
    Group = case_when(
      short_prefix == "ASO" ~ "CTL",
      short_prefix == "MCI" ~ "MCI",
      short_prefix == "AD"  ~ "AD",
      TRUE                  ~ "DROP_ASM" 
    )
  ) %>%
  # Filter out ASM samples entirely from the clinical layout
  filter(Group != "DROP_ASM")

metadata_rna$Group <- factor(metadata_rna$Group, levels = c("CTL", "MCI", "AD"))

# ==========================================
# STEP 3: Load and Synchronize Count Matrix
# ==========================================
cat("Loading local compressed counts matrix file...\n")
matrix_path <- "GSE270454/GSE270454_RNAseq-combined-counts-matrix.csv.gz"

if(!file.exists(matrix_path)) {
  stop("Missing file: Please ensure 'GSE270454/GSE270454_RNAseq-combined-counts-matrix.csv.gz' is in your directory.")
}

counts_all <- read.csv(matrix_path, row.names = 1, check.names = FALSE)

# Slice out only the 34 active target columns matching your filtered metadata rows
rna_sample_ids <- metadata_rna$count_col_id
counts_rna     <- counts_all[, rna_sample_ids]

rownames(metadata_rna) <- metadata_rna$count_col_id
metadata_rna           <- metadata_rna[colnames(counts_rna), ]

# ==========================================
# STEP 4: Run the Statistical Model
# ==========================================
cat("Initializing and running DESeq2 model...\n")
dds  <- DESeqDataSetFromMatrix(countData = counts_rna, colData = metadata_rna, design = ~ Group)
keep <- rowSums(counts(dds) >= 10) >= 10 
dds  <- dds[keep,]
dds  <- DESeq(dds)

# ==========================================
# FILE 1: Log2-Transformed Expression Matrix Generation
# ==========================================
cat("Generating File 1: Log2-Transformed Expression Matrix...\n")
expression_matrix <- counts(dds, normalized = TRUE)

# Apply log2 conversion with a pseudo-count of 1 to preserve zeros safely
log2_expression_matrix <- log2(expression_matrix + 1)

# Map columns directly to official GSM accession numbers (e.g., GSM8343242)
colnames(log2_expression_matrix) <- metadata_rna$geo_accession[match(colnames(log2_expression_matrix), metadata_rna$count_col_id)]

write.csv(log2_expression_matrix, "RNAseq_Normalized_Expression_Matrix.csv", row.names = TRUE)

# ==========================================
# FILE 2 & 3: Custom DEG Processing
# ==========================================
cat("Extracting statistical model outputs...\n")

res_MCI_vs_CTL <- as.data.frame(results(dds, contrast = c("Group", "MCI", "CTL")))
res_AD_vs_CTL  <- as.data.frame(results(dds, contrast = c("Group", "AD", "CTL")))
res_MCI_vs_AD  <- as.data.frame(results(dds, contrast = c("Group", "MCI", "AD"))) 

# Standardize layout to your requested 6 columns
standardize_columns <- function(df, group_label) {
  df_clean <- df %>% rownames_to_column(var = "gene_symbol")
  
  colnames(df_clean)[grepl("log2FoldChange", colnames(df_clean))] <- "logFC"
  colnames(df_clean)[grepl("^stat$", colnames(df_clean))]         <- "t_value"
  colnames(df_clean)[grepl("^pvalue$", colnames(df_clean))]       <- "P_value"
  colnames(df_clean)[grepl("^padj$", colnames(df_clean))]         <- "adj_p_value"
  
  df_clean %>%
    mutate(Group = group_label) %>%
    select(gene_symbol, logFC, P_value, adj_p_value, t_value, Group)
}

deg_MCI_vs_CTL <- standardize_columns(res_MCI_vs_CTL, "MCI_vs_CTL")
deg_AD_vs_CTL  <- standardize_columns(res_AD_vs_CTL, "AD_vs_CTL")
deg_MCI_vs_AD  <- standardize_columns(res_MCI_vs_AD, "MCI_vs_AD")

# Bind all rows completely into the master structure
deg_combined <- rbind(deg_MCI_vs_CTL, deg_AD_vs_CTL, deg_MCI_vs_AD)

# Write File 2: Completely Unfiltered Master File
cat("Writing File 2: 'RNAseq_Combined_DEG_Analysis.csv' (Unfiltered)...\n")
write.csv(deg_combined, "RNAseq_Combined_DEG_Analysis.csv", row.names = FALSE)

# Filter File 3: Enforce strict Microarray settings (p < 0.05, |logFC| > 0.2)
sig_deg <- deg_combined %>%
  filter(P_value < 0.05) %>%
  filter(abs(logFC) > 0.2)

# Write File 3: Screened Validating Signatures
cat("Writing File 3: 'sigDEG.csv' (Strict Filtered Targets)...\n")
write.csv(sig_deg, "sigDEG.csv", row.names = FALSE)

# ==========================================
# DATA INTEGRITY VERIFICATION SUMMARY
# ==========================================
cat("\n============================================\n")
cat("          PIPELINE EXECUTION SUCCESS        \n")
cat("============================================\n")
cat("Processed Samples        :", ncol(counts_rna), "\n") # Will print 34
cat("Matrix/Metadata Sync     :", all(colnames(log2_expression_matrix) == metadata_rna$geo_accession), "\n")
cat("Total Complete DEG Rows  :", nrow(deg_combined), "\n")
cat("Significant Validated DEGs:", nrow(sig_deg), "rows passed limits.\n\n")
print(table(sig_deg$Group))
cat("\nAll 3 files updated successfully with ASO controls only!\n")
cat("============================================\n")
