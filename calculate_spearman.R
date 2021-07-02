## To evaluate the correlation between each cdks_ccns, I calculate the Spearman correlation coefficient.
## Since gene expression does not follow normal distribution, The spearsman is more appropriate than Pearson.


# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Load package
library(tidyverse)
library(magrittr)
library(doParallel)

# Define cdks_ccns_list and cancer_list
basedir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data"
savedir <- "/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer"
cdks_ccns <- read_tsv(file.path(savedir, "cdks_ccns.txt"))[[1]]
setwd(basedir)
cancers <- list.files()
cancers %<>%
  .[-which(. %in% c("COAD_READ_icluding_normal_tissue_nature_2012", "COAD_READ_nonreference",
    "GBM_nature_2008", "GBM_nonreference"))]

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

cal_corrcoeff <- function(target_gene) {
  target_gene_tpm <- tpm %>%
    filter(., rownames(.) == target_gene)
  df_tpm_cor <- tpm
  # Caluculate spearman correlation coefficient
  cordata_tpm <- apply(df_tpm_cor, 1, function(x) {
    cor.test(as.numeric(target_gene_tpm), as.numeric(x), method = "spearman")
  })

  # Join the rho_value and p_value
  rho_tpm <- sapply(cordata_tpm, "[[", "estimate")
  pvalue_tpm <- sapply(cordata_tpm, "[[", "p.value")
  name_rho_tpm <- gsub(names(rho_tpm), pattern = ".rho", replacement = "")
  dif_name <- identical(name_rho_tpm, names(pvalue_tpm))
  if (dif_name == TRUE) {
    cor_tpm <- cbind(rho_tpm, pvalue_tpm)
  }

  # Remove '.rho' from rowname
  rownames(cor_tpm) <- gsub(rownames(cor_tpm), pattern = ".rho", replacement = "")

  # Save data
  setwd(file.path(savedir, "Datasets"))
  make_folder("Cor_Datasets")
  setwd("Cor_Datasets")
  make_folder(toupper(target_cancer))
  setwd(toupper(target_cancer))
  cor_tpm_save <- as_tibble(rownames_to_column(as.data.frame(cor_tpm),
    "gene_symbol"))
  write_tsv(cor_tpm_save, paste0(tolower(target_gene), "_cor_tpm.txt"))
}

# Ignores the warning message: Cannot compute exact p-value with ties
options(warn = -1)

for (i in 1:length(cancers)) {
  # Define target_cancer
  target_cancer <- cancers[i]

  # Load tpm data
  setwd(file.path(basedir, toupper(target_cancer)))
  txt_tpm <- list.files()[grep("_tpm\\.txt$", list.files())]
  tpm <- read_tsv(txt_tpm) %>%
    .[, -which(colnames(.) == "ensembl_id")] %>%
    as.data.frame(.)
  rownames(tpm) <- tpm[["gene_symbol"]]
  tpm %<>%
    .[, -which(colnames(.) == "gene_symbol")]

  # Remove low expression gene
  param_thres <- 0
  obj <- as.logical(rowSums(tpm) > param_thres)
  tpm %<>%
    .[obj, ]

  # Reconstruct cdks_ccns because low expression genes are removed
  cdks_ccns %<>%
    intersect(., rownames(tpm))

  # Execute parallel processing
  cores <- getOption("mc.cores", detectCores())
  cl <- makeCluster(cores, type = "FORK")
  registerDoParallel(cl)
  parLapply(cl, cdks_ccns, cal_corrcoeff)
  stopCluster(cl)
}

traceback()

