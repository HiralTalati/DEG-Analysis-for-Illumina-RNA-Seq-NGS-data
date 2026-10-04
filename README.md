**RNA-Seq Differential Expression Analysis Pipeline (DESeq2)**

This repository contains an end-to-end transcriptomic pipeline written in R for analyzing Next-Generation Sequencing (NGS) data generated on the Illumina platform. The pipeline processes raw gene counts using DESeq2, performs sample harmonization against clinical metadata from NCBI GEO, and runs multi-group statistical comparisons for staging disease progression (Control vs. MCI vs. AD).

📌 Features
• Automated Environment Setup: Checks, installs, and mounts missing CRAN and Bioconductor packages (DESeq2, tidyverse, GEOquery).

• Robust Hybrid Metadata Fetching: Retrieves sample characteristics using an online API, with an automatic fall-back mechanism to scan local file directories if connection errors occur.

• Strict Age/Condition Harmonization: Standardises raw clinical titles into precise diagnostic cohorts. Specifically drops middle-aged/older adults (ASM) to restrict the baseline control group exclusively to older adults (ASO → CTL), alongside MCI and AD patient groups.

• Low-Count Row Filtering: Filters out low-abundance genes, retaining only features with ≥ 10 counts in ≥ 10 samples to boost statistical testing power.

• Multi-Group Contrast Modeling: Executes three pairwise group comparisons to evaluate transition dynamics (MCI vs. CTL, AD vs. CTL, and MCI vs. AD).

• Dual Statistical Reshaping: Formats output matrices to match conventional downstream analysis layouts, standardizing variable headers to logFC, t_value, P_value, and adj_p_value.

🛠️ Tech Stack & Dependencies

The pipeline is fully compatible with R (v4.0+) environments.

Required Libraries

# Checked and installed automatically at runtime
BiocManager, GEOquery, tidyverse, DESeq2

📂 Required Workspace Structure
The script operates on the GSE270454 dataset index. Before launching the execution cycle, ensure your workspace folder is set up as follows:

├── GSE270454/

│   ├── GSE270454_series_matrix.txt.gz           # Optional: local metadata fallback


│   └── GSE270454_RNAseq-combined-counts-matrix.csv.gz # REQUIRED: Compressed raw NGS counts

├── deseq2_pipeline.R                            # This script file

└── README.md


🚀 Pipeline Workflow
1. Synchronized Matrix Slicing

The raw count matrix is loaded and matched against the filtered metadata rows. The script filters out columns corresponding to the dropped ASM cohort, extracting only the 34 active target columns that match the strict aging parameters.

2. Running the Model
   
The generalized linear model (GLM) is fitted using a single-factor design:

design = ~ Group

Following sample filtering and library size normalization, differential expression is calculated using negative binomial distribution parameters via the master wrapper DESeq(dds).

📄 License
This analysis framework is open-source and distributed under the MIT License.
