# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Load package
library(tidyverse)
library(recount3)
library(recount)
library(SummarizedExperiment)
library(doParallel)

# Define parameters
basedir <- "/Volumes/G_DRIVEmobile/recount3_Rdata"
human_projects <- available_projects()
human_projects_gtex <- subset(human_projects, file_source == "gtex" & project_type ==
  "data_sources")
human_projects_tcga <- subset(human_projects, file_source == "tcga" & project_type ==
  "data_sources")

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_rse <- function(target, file_source) {
  # Define parameters
  setwd(basedir)
  make_folder(file_source)
  basedir <- file.path(basedir, file_source)
  setwd(basedir)
  make_folder(toupper(target))
  save_dir <- file.path(basedir, toupper(target))

  # Load rse_gene
  proj_info <- human_projects %>%
    subset(., file_source == tolower(file_source) & project == toupper(target))
  rse_gene <- create_rse(proj_info)

  # Scale the counts to a common denominator (default = 40M)
  # Default option: AUC
  assay(rse_gene, "counts") <- transform_counts(rse_gene)

  # Calculate TPM and put data in rse_object
  assays(rse_gene)$TPM <- recount::getTPM(rse_gene, length_var = "score")

  # Make data_frame of counts and TPM
  df_counts <- as.data.frame(assay(rse_gene, "counts")) %>%
    add_column(., ensembl_id = rownames(.), .before = colnames(.)[1])
  df_tpm <- as.data.frame(assay(rse_gene, "TPM")) %>%
    add_column(., ensembl_id = rownames(.), .before = colnames(.)[1])
  df_meta <- as.data.frame(colData(rse_gene))

  # Replace the rownames to gene symbol
  rownames(df_counts) <- make.unique(rowData(rse_gene)$gene_name)
  df_counts <- add_column(df_counts, gene_symbol = rownames(df_counts),
    .after = colnames(df_counts)[1])
  rownames(df_tpm) <- make.unique(rowData(rse_gene)$gene_name)
  df_tpm <- add_column(df_tpm, gene_symbol = rownames(df_tpm), .after = colnames(df_tpm)[1])

  # Save files
  save_counts <- file.path(save_dir, paste0(target, "_counts.txt"))
  save_tpm <- file.path(save_dir, paste0(target, "_tpm.txt"))
  save_meta <- file.path(save_dir, paste0(target, "_meta.txt"))
  write_tsv(df_counts, save_counts)
  write_tsv(df_tpm, save_tpm)
  write_tsv(df_meta, save_meta)
  save_rse_gene <- file.path(save_dir, paste0(target, ".rda"))
  save(rse_gene, file = save_rse_gene)
}

# Make project lists
projects_gtex <- tolower(human_projects_gtex$project)
projects_tcga <- tolower(human_projects_tcga$project)

# Execute parallel processing
cores <- getOption("mc.cores", detectCores())
cl <- makeCluster(cores, type = "FORK")
registerDoParallel(cl)
parLapply(cl, projects_gtex, make_rse, "GTEx")
# Stop cluster
stopCluster(cl)

# Execute parallel processing
cores <- getOption("mc.cores", detectCores())
cl <- makeCluster(cores, type = "FORK")
registerDoParallel(cl)
parLapply(cl, projects_tcga, make_rse, "TCGA")
# Stop cluster
stopCluster(cl)

traceback()

