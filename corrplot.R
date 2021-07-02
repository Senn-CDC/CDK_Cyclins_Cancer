## Based on the Spearman correlation coefficient, correlate heatmaps are plotted for each cancer

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Load package
library(tidyverse)
library(magrittr)
library(corrplot)

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

# Make Cor_Heatmap
for (i in 1:length(cancers)) {
  target_cancer <- cancers[i]
  df_corrplot <- tibble(gene_symbol = cdks_ccns)
  # Load corralation data
  for (j in 1:length(cdks_ccns)) {
    target_cdks_ccns <- cdks_ccns[j]
    basedir <- "/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/Datasets/Cor_Datasets"
    target_dir <- file.path(basedir, toupper(target_cancer))
    setwd(target_dir)
    if (target_cdks_ccns %in% toupper(gsub("_cor_tpm\\.txt", "", list.files()))) {
      txt_tpm <- list.files()[grep("_tpm\\.txt$", list.files())]
      txt_corr <- read_tsv(file.path(target_dir, paste0(target_cdks_ccns,
        "_cor_tpm.txt")))
      # Select CDKs and CCNs data
      txt_corr %<>%
        filter_at(., vars(all_of("gene_symbol")), all_vars(. %in%
          cdks_ccns)) %>%
        mutate(., gene_symbol = factor(gene_symbol, levels = cdks_ccns)) %>%
        arrange(gene_symbol)
      txt_corr %<>%
        select(., c("gene_symbol", "rho_tpm"))
      df_corrplot <- inner_join(df_corrplot, txt_corr, by = "gene_symbol")
      df_corrplot %<>%
        rename_at(., vars(all_of("rho_tpm")), ~target_cdks_ccns)
    } else {
      next
    }
  }

  df_corrplot %<>%
    as.data.frame(.)
  rownames(df_corrplot) <- df_corrplot[["gene_symbol"]]
  df_corrplot %<>%
    .[, -which(colnames(.) == "gene_symbol")]

  # Define output parameter
  output_height <- 15
  output_width <- 15

  # Save data
  setwd(savedir)
  make_folder("Cor_Heatmap")
  setwd("Cor_Heatmap")
  make_folder("CDKs_CCNs")
  setwd("CDKs_CCNs")
  pdf(file = paste0(tolower(target_cancer), "_spearman_correlation_coefficient.pdf"),
    height = output_height, width = output_width)
  corrplot(as.matrix(df_corrplot), method = "circle", title = paste0(target_cancer,
    "\nSpearman Correlation Coefficient"), type = "upper", tl.col = "black",
    tl.cex = 0.6, tl.srt = 45, mar = c(0, 0, 7, 0), cl.cex = 1.2)
  dev.off()
}

traceback()

