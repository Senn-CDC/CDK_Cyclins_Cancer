## To evaluate gene expression levels, I make two types of heatmap, ratio heatmap and clustred heatmap
## The ratio heatmap: To evaluate gene expression levels for each cancer, I adjusted the width of t
## he heatmap and line them up from high expression to low expression.
## So sample's order are not corresponding to each CDKs and CCNs.
## Clustered_heatmap: To evaluate the relationship for each CDKs and CCNs, I perfumed clustering each sample. 
## So sample's order are corresponding to each CDKs and CCNs

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Load package
library(tidyverse)
library(magrittr)
library(ComplexHeatmap)
library(circlize)
library(dendextend)

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

# Load data and select cdks_ccns samples
df_ht_tpm <- read_tsv("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/Datasets/tcga_cdks_ccns_tpm.txt")

# Define cancer_list
cancers <- unique(gsub("\\d", "", colnames(df_ht_tpm))) %>%
  .[-which(. == "gene_symbol")]

# Make Heatmap
# -------------------------------------------------------------------------
# Define Save directry
save_predir <- file.path("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer")
setwd(save_predir)
make_folder("Expression_Heatmap")
save_dir <- file.path(save_predir, "Expression_Heatmap")
setwd(save_dir)
make_folder("tpm")
setwd(file.path(save_dir, "tpm"))

# Ignores the warning message: Heatmap/annotation names are
# duplicated:Heatmap/annotation names are duplicated: gene expression
options(warn = -1)

# Define legend_param
legend_title <- "gene expression (TPM)"
bottom <- 0
top <- 50
middle <- top/2
col_fun <- colorRamp2(c(bottom, middle, top), c("blue", "white", "red"))
heatmap_legend_param <- list(title_position = "leftcenter-rot", legend_height = unit(18,
  "cm"), grid_width = unit(1, "cm"), title_gp = gpar(fontsize = 20),
  labels_gp = gpar(fontsize = 20, font = 10))

# Define cdks_ccns_list
cdks_ccns <- read_tsv("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/cdks_ccns.txt")[[1]]

# Reconstract dataframe
df_cdks_ccns_ht_tpm <- df_ht_tpm %>%
  filter_at(., vars(all_of("gene_symbol")), all_vars(. %in% cdks_ccns)) %>%
  mutate(., gene_symbol = factor(gene_symbol, levels = cdks_ccns)) %>%
  arrange(gene_symbol)
df_cdks_ccns_ht_tpm <- as.data.frame(df_cdks_ccns_ht_tpm)
rownames(df_cdks_ccns_ht_tpm) <- df_cdks_ccns_ht_tpm[["gene_symbol"]]
df_cdks_ccns_ht_tpm %<>%
  .[, -which(colnames(.) == "gene_symbol")]

# Setting the global parameter for heatmap title
ht_opt$TITLE_PADDING <- unit(c(8.5, 8.5), "points")

# Define output parameter
output_height <- 25
output_width <- 45

ht_tpm <- NULL
for (i in 1:length(cancers)) {
  target_cancer <- cancers[i]
  df_target_cancer_ht_tpm <- select(df_cdks_ccns_ht_tpm, starts_with(target_cancer))
  # cluster_cols: Calculate distance between each sample, and put the
  # result into df_target_cancer_ht_tpm_dist_column
  c_column <- dist(t(df_target_cancer_ht_tpm), method = "euclidean")
  df_target_cancer_ht_tpm_dist_column <- as.dist(c_column)
  # Execute clustering
  column_dend <- hclust(df_target_cancer_ht_tpm_dist_column, method = "average")
  # Define heatmap_param
  column_title <- paste0(target_cancer, "\n", "(", length(df_target_cancer_ht_tpm),
    ")")
  column_order <- colnames(df_target_cancer_ht_tpm)
  column_title_gp <- gpar(fontsize = 5, fontface = "bold", fill = i,
    col = "white", border = "white")
  row_title <- " "
  pre_ht_tpm <- Heatmap(as.matrix(df_target_cancer_ht_tpm), name = legend_title,
    col = col_fun, cluster_columns = column_dend, column_title = column_title,
    column_title_gp = column_title_gp, show_column_dend = FALSE, show_column_names = FALSE,
    row_split = 1:length(cdks_ccns), row_order = cdks_ccns, row_title = row_title,
    row_names_gp = gpar(fontsize = 10), row_names_side = "left", show_row_dend = FALSE,
    show_row_names = TRUE, height = unit(50, "cm"), heatmap_legend_param = heatmap_legend_param,
    use_raster = TRUE)
  ht_tpm <- ht_tpm + pre_ht_tpm
  pdf(paste0(top, "_CDKs_CCNs", "_Clustered_Heatmap.pdf"), height = output_height,
    width = output_width)
  draw(ht_tpm)
  dev.off()
}

# Define cdks_ccns_list
cdks_ccns <- read_tsv("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/cdks.txt")[[1]]

# Reconstract dataframe
df_cdks_ccns_ht_tpm <- df_ht_tpm %>%
  filter_at(., vars(all_of("gene_symbol")), all_vars(. %in% cdks_ccns)) %>%
  mutate(., gene_symbol = factor(gene_symbol, levels = cdks_ccns)) %>%
  arrange(gene_symbol)
df_cdks_ccns_ht_tpm <- as.data.frame(df_cdks_ccns_ht_tpm)
rownames(df_cdks_ccns_ht_tpm) <- df_cdks_ccns_ht_tpm[["gene_symbol"]]
df_cdks_ccns_ht_tpm %<>%
  .[, -which(colnames(.) == "gene_symbol")]


# Setting the global parameter for heatmap title
ht_opt$TITLE_PADDING <- unit(c(8.5, 8.5), "points")

# Define output parameter
output_height <- 15
output_width <- 60

ht_tpm <- NULL
for (i in 1:length(cancers)) {
  target_cancer <- cancers[i]
  df_target_cancer_ht_tpm <- select(df_cdks_ccns_ht_tpm, starts_with(target_cancer))
  # cluster_cols: Calculate distance between each sample, and put the
  # result into df_target_cancer_ht_tpm_dist_column
  c_column <- dist(t(df_target_cancer_ht_tpm), method = "euclidean")
  df_target_cancer_ht_tpm_dist_column <- as.dist(c_column)
  # Execute clustering
  column_dend <- hclust(df_target_cancer_ht_tpm_dist_column, method = "average")
  # Define heatmap_param
  column_title <- paste0(target_cancer, "\n", "(", length(df_target_cancer_ht_tpm),
    ")")
  column_order <- colnames(df_target_cancer_ht_tpm)
  column_title_gp <- gpar(fontsize = 5, fontface = "bold", fill = i,
    col = "white", border = "white")
  row_title <- " "
  pre_ht_tpm <- Heatmap(as.matrix(df_target_cancer_ht_tpm), name = legend_title,
    col = col_fun, cluster_columns = column_dend, column_title = column_title,
    column_title_gp = column_title_gp, show_column_dend = FALSE, show_column_names = FALSE,
    row_split = 1:length(cdks_ccns), row_order = cdks_ccns, row_title = row_title,
    row_names_gp = gpar(fontsize = 12), row_names_side = "left", show_row_dend = FALSE,
    show_row_names = TRUE, height = unit(25, "cm"), heatmap_legend_param = heatmap_legend_param,
    use_raster = TRUE)
  ht_tpm <- ht_tpm + pre_ht_tpm
  pdf(paste0(top, "_CDKs", "_Clustered_Heatmap.pdf"), height = output_height,
    width = output_width)
  draw(ht_tpm)
  dev.off()
}

# Define cdks_ccns_list
cdks_ccns <- read_tsv("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/cdks.txt")[[1]]

# Reconstract dataframe
df_cdks_ccns_ht_tpm <- df_ht_tpm %>%
  filter_at(., vars(all_of("gene_symbol")), all_vars(. %in% cdks_ccns)) %>%
  mutate(., gene_symbol = factor(gene_symbol, levels = cdks_ccns)) %>%
  arrange(gene_symbol)
df_cdks_ccns_ht_tpm <- as.data.frame(df_cdks_ccns_ht_tpm)
rownames(df_cdks_ccns_ht_tpm) <- df_cdks_ccns_ht_tpm[["gene_symbol"]]
df_cdks_ccns_ht_tpm %<>%
  .[, -which(colnames(.) == "gene_symbol")]


# Setting the global parameter for heatmap title
ht_opt$TITLE_PADDING <- unit(c(8.5, 8.5), "points")

# Define output parameter
output_height <- 15
output_width <- 60

ht_tpm <- NULL
for (i in 1:length(cancers)) {
  target_cancer <- cancers[i]
  df_target_cancer_ht_tpm <- select(df_cdks_ccns_ht_tpm, starts_with(target_cancer))
  # cluster_cols: Calculate distance between each sample, and put the
  # result into df_target_cancer_ht_tpm_dist_column
  c_column <- dist(t(df_target_cancer_ht_tpm), method = "euclidean")
  df_target_cancer_ht_tpm_dist_column <- as.dist(c_column)
  # Execute clustering
  column_dend <- hclust(df_target_cancer_ht_tpm_dist_column, method = "average")
  # Define heatmap_param
  column_title <- paste0(target_cancer, "\n", "(", length(df_target_cancer_ht_tpm),
    ")")
  column_order <- colnames(df_target_cancer_ht_tpm)
  column_title_gp <- gpar(fontsize = 5, fontface = "bold", fill = i,
    col = "white", border = "white")
  row_title <- " "
  pre_ht_tpm <- Heatmap(as.matrix(df_target_cancer_ht_tpm), name = legend_title,
    col = col_fun, cluster_columns = column_dend, column_title = column_title,
    column_title_gp = column_title_gp, show_column_dend = FALSE, show_column_names = FALSE,
    row_split = 1:length(cdks_ccns), row_order = cdks_ccns, row_title = row_title,
    row_names_gp = gpar(fontsize = 11), row_names_side = "left", show_row_dend = FALSE,
    show_row_names = TRUE, height = unit(25, "cm"), heatmap_legend_param = heatmap_legend_param,
    use_raster = TRUE)
  ht_tpm <- ht_tpm + pre_ht_tpm
  pdf(paste0(top, "_CCNs", "_Clustered_Heatmap.pdf"), height = output_height,
    width = output_width)
  draw(ht_tpm)
  dev.off()
}

# Define legend_param
legend_title <- "gene expression (TPM)"
bottom <- 0
top <- 100
middle <- top/2
col_fun <- colorRamp2(c(bottom, middle, top), c("blue", "white", "red"))
heatmap_legend_param <- list(title_position = "leftcenter-rot", legend_height = unit(18,
  "cm"), grid_width = unit(1, "cm"), title_gp = gpar(fontsize = 20),
  labels_gp = gpar(fontsize = 20, font = 10))

# Define cdks_ccns_list
cdks_ccns <- read_tsv("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/cdks_ccns.txt")[[1]]

# Reconstract dataframe
df_cdks_ccns_ht_tpm <- df_ht_tpm %>%
  filter_at(., vars(all_of("gene_symbol")), all_vars(. %in% cdks_ccns)) %>%
  mutate(., gene_symbol = factor(gene_symbol, levels = cdks_ccns)) %>%
  arrange(gene_symbol)
df_cdks_ccns_ht_tpm <- as.data.frame(df_cdks_ccns_ht_tpm)
rownames(df_cdks_ccns_ht_tpm) <- df_cdks_ccns_ht_tpm[["gene_symbol"]]
df_cdks_ccns_ht_tpm %<>%
  .[, -which(colnames(.) == "gene_symbol")]

# Setting the global parameter for heatmap title
ht_opt$TITLE_PADDING <- unit(c(8.5, 8.5), "points")

# Define output parameter
output_height <- 25
output_width <- 45

ht_tpm <- NULL
for (i in 1:length(cancers)) {
  target_cancer <- cancers[i]
  df_target_cancer_ht_tpm <- select(df_cdks_ccns_ht_tpm, starts_with(target_cancer))
  # cluster_cols: Calculate distance between each sample, and put the
  # result into df_target_cancer_ht_tpm_dist_column
  c_column <- dist(t(df_target_cancer_ht_tpm), method = "euclidean")
  df_target_cancer_ht_tpm_dist_column <- as.dist(c_column)
  # Execute clustering
  column_dend <- hclust(df_target_cancer_ht_tpm_dist_column, method = "average")
  # Define heatmap_param
  column_title <- paste0(target_cancer, "\n", "(", length(df_target_cancer_ht_tpm),
    ")")
  column_order <- colnames(df_target_cancer_ht_tpm)
  column_title_gp <- gpar(fontsize = 5, fontface = "bold", fill = i,
    col = "white", border = "white")
  row_title <- " "
  pre_ht_tpm <- Heatmap(as.matrix(df_target_cancer_ht_tpm), name = legend_title,
    col = col_fun, cluster_columns = column_dend, column_title = column_title,
    column_title_gp = column_title_gp, show_column_dend = FALSE, show_column_names = FALSE,
    row_split = 1:length(cdks_ccns), row_order = cdks_ccns, row_title = row_title,
    row_names_gp = gpar(fontsize = 10), row_names_side = "left", show_row_dend = FALSE,
    show_row_names = TRUE, height = unit(50, "cm"), heatmap_legend_param = heatmap_legend_param,
    use_raster = TRUE)
  ht_tpm <- ht_tpm + pre_ht_tpm
  pdf(paste0(top, "_CDKs_CCNs", "_Clustered_Heatmap.pdf"), height = output_height,
    width = output_width)
  draw(ht_tpm)
  dev.off()
}

# Define cdks_ccns_list
cdks_ccns <- read_tsv("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/cdks.txt")[[1]]

# Reconstract dataframe
df_cdks_ccns_ht_tpm <- df_ht_tpm %>%
  filter_at(., vars(all_of("gene_symbol")), all_vars(. %in% cdks_ccns)) %>%
  mutate(., gene_symbol = factor(gene_symbol, levels = cdks_ccns)) %>%
  arrange(gene_symbol)
df_cdks_ccns_ht_tpm <- as.data.frame(df_cdks_ccns_ht_tpm)
rownames(df_cdks_ccns_ht_tpm) <- df_cdks_ccns_ht_tpm[["gene_symbol"]]
df_cdks_ccns_ht_tpm %<>%
  .[, -which(colnames(.) == "gene_symbol")]


# Setting the global parameter for heatmap title
ht_opt$TITLE_PADDING <- unit(c(8.5, 8.5), "points")

# Define output parameter
output_height <- 15
output_width <- 60

ht_tpm <- NULL
for (i in 1:length(cancers)) {
  target_cancer <- cancers[i]
  df_target_cancer_ht_tpm <- select(df_cdks_ccns_ht_tpm, starts_with(target_cancer))
  # cluster_cols: Calculate distance between each sample, and put the
  # result into df_target_cancer_ht_tpm_dist_column
  c_column <- dist(t(df_target_cancer_ht_tpm), method = "euclidean")
  df_target_cancer_ht_tpm_dist_column <- as.dist(c_column)
  # Execute clustering
  column_dend <- hclust(df_target_cancer_ht_tpm_dist_column, method = "average")
  # Define heatmap_param
  column_title <- paste0(target_cancer, "\n", "(", length(df_target_cancer_ht_tpm),
    ")")
  column_order <- colnames(df_target_cancer_ht_tpm)
  column_title_gp <- gpar(fontsize = 5, fontface = "bold", fill = i,
    col = "white", border = "white")
  row_title <- " "
  pre_ht_tpm <- Heatmap(as.matrix(df_target_cancer_ht_tpm), name = legend_title,
    col = col_fun, cluster_columns = column_dend, column_title = column_title,
    column_title_gp = column_title_gp, show_column_dend = FALSE, show_column_names = FALSE,
    row_split = 1:length(cdks_ccns), row_order = cdks_ccns, row_title = row_title,
    row_names_gp = gpar(fontsize = 12), row_names_side = "left", show_row_dend = FALSE,
    show_row_names = TRUE, height = unit(25, "cm"), heatmap_legend_param = heatmap_legend_param,
    use_raster = TRUE)
  ht_tpm <- ht_tpm + pre_ht_tpm
  pdf(paste0(top, "_CDKs", "_Clustered_Heatmap.pdf"), height = output_height,
    width = output_width)
  draw(ht_tpm)
  dev.off()
}

# Define cdks_ccns_list
cdks_ccns <- read_tsv("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/cdks.txt")[[1]]

# Reconstract dataframe
df_cdks_ccns_ht_tpm <- df_ht_tpm %>%
  filter_at(., vars(all_of("gene_symbol")), all_vars(. %in% cdks_ccns)) %>%
  mutate(., gene_symbol = factor(gene_symbol, levels = cdks_ccns)) %>%
  arrange(gene_symbol)
df_cdks_ccns_ht_tpm <- as.data.frame(df_cdks_ccns_ht_tpm)
rownames(df_cdks_ccns_ht_tpm) <- df_cdks_ccns_ht_tpm[["gene_symbol"]]
df_cdks_ccns_ht_tpm %<>%
  .[, -which(colnames(.) == "gene_symbol")]


# Setting the global parameter for heatmap title
ht_opt$TITLE_PADDING <- unit(c(8.5, 8.5), "points")

# Define output parameter
output_height <- 15
output_width <- 60

ht_tpm <- NULL
for (i in 1:length(cancers)) {
  target_cancer <- cancers[i]
  df_target_cancer_ht_tpm <- select(df_cdks_ccns_ht_tpm, starts_with(target_cancer))
  # cluster_cols: Calculate distance between each sample, and put the
  # result into df_target_cancer_ht_tpm_dist_column
  c_column <- dist(t(df_target_cancer_ht_tpm), method = "euclidean")
  df_target_cancer_ht_tpm_dist_column <- as.dist(c_column)
  # Execute clustering
  column_dend <- hclust(df_target_cancer_ht_tpm_dist_column, method = "average")
  # Define heatmap_param
  column_title <- paste0(target_cancer, "\n", "(", length(df_target_cancer_ht_tpm),
    ")")
  column_order <- colnames(df_target_cancer_ht_tpm)
  column_title_gp <- gpar(fontsize = 5, fontface = "bold", fill = i,
    col = "white", border = "white")
  row_title <- " "
  pre_ht_tpm <- Heatmap(as.matrix(df_target_cancer_ht_tpm), name = legend_title,
    col = col_fun, cluster_columns = column_dend, column_title = column_title,
    column_title_gp = column_title_gp, show_column_dend = FALSE, show_column_names = FALSE,
    row_split = 1:length(cdks_ccns), row_order = cdks_ccns, row_title = row_title,
    row_names_gp = gpar(fontsize = 11), row_names_side = "left", show_row_dend = FALSE,
    show_row_names = TRUE, height = unit(25, "cm"), heatmap_legend_param = heatmap_legend_param,
    use_raster = TRUE)
  ht_tpm <- ht_tpm + pre_ht_tpm
  pdf(paste0(top, "_CCNs", "_Clustered_Heatmap.pdf"), height = output_height,
    width = output_width)
  draw(ht_tpm)
  dev.off()
}

# Define legend_param
legend_title <- "gene expression (TPM)"
bottom <- 0
top <- 150
middle <- top/2
col_fun <- colorRamp2(c(bottom, middle, top), c("blue", "white", "red"))
heatmap_legend_param <- list(title_position = "leftcenter-rot", legend_height = unit(18,
  "cm"), grid_width = unit(1, "cm"), title_gp = gpar(fontsize = 20),
  labels_gp = gpar(fontsize = 20, font = 10))

# Define cdks_ccns_list
cdks_ccns <- read_tsv("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/cdks_ccns.txt")[[1]]

# Reconstract dataframe
df_cdks_ccns_ht_tpm <- df_ht_tpm %>%
  filter_at(., vars(all_of("gene_symbol")), all_vars(. %in% cdks_ccns)) %>%
  mutate(., gene_symbol = factor(gene_symbol, levels = cdks_ccns)) %>%
  arrange(gene_symbol)
df_cdks_ccns_ht_tpm <- as.data.frame(df_cdks_ccns_ht_tpm)
rownames(df_cdks_ccns_ht_tpm) <- df_cdks_ccns_ht_tpm[["gene_symbol"]]
df_cdks_ccns_ht_tpm %<>%
  .[, -which(colnames(.) == "gene_symbol")]

# Setting the global parameter for heatmap title
ht_opt$TITLE_PADDING <- unit(c(8.5, 8.5), "points")

# Define output parameter
output_height <- 25
output_width <- 45

ht_tpm <- NULL
for (i in 1:length(cancers)) {
  target_cancer <- cancers[i]
  df_target_cancer_ht_tpm <- select(df_cdks_ccns_ht_tpm, starts_with(target_cancer))
  # cluster_cols: Calculate distance between each sample, and put the
  # result into df_target_cancer_ht_tpm_dist_column
  c_column <- dist(t(df_target_cancer_ht_tpm), method = "euclidean")
  df_target_cancer_ht_tpm_dist_column <- as.dist(c_column)
  # Execute clustering
  column_dend <- hclust(df_target_cancer_ht_tpm_dist_column, method = "average")
  # Define heatmap_param
  column_title <- paste0(target_cancer, "\n", "(", length(df_target_cancer_ht_tpm),
    ")")
  column_order <- colnames(df_target_cancer_ht_tpm)
  column_title_gp <- gpar(fontsize = 5, fontface = "bold", fill = i,
    col = "white", border = "white")
  row_title <- " "
  pre_ht_tpm <- Heatmap(as.matrix(df_target_cancer_ht_tpm), name = legend_title,
    col = col_fun, cluster_columns = column_dend, column_title = column_title,
    column_title_gp = column_title_gp, show_column_dend = FALSE, show_column_names = FALSE,
    row_split = 1:length(cdks_ccns), row_order = cdks_ccns, row_title = row_title,
    row_names_gp = gpar(fontsize = 10), row_names_side = "left", show_row_dend = FALSE,
    show_row_names = TRUE, height = unit(50, "cm"), heatmap_legend_param = heatmap_legend_param,
    use_raster = TRUE)
  ht_tpm <- ht_tpm + pre_ht_tpm
  pdf(paste0(top, "_CDKs_CCNs", "_Clustered_Heatmap.pdf"), height = output_height,
    width = output_width)
  draw(ht_tpm)
  dev.off()
}

# Define cdks_ccns_list
cdks_ccns <- read_tsv("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/cdks.txt")[[1]]

# Reconstract dataframe
df_cdks_ccns_ht_tpm <- df_ht_tpm %>%
  filter_at(., vars(all_of("gene_symbol")), all_vars(. %in% cdks_ccns)) %>%
  mutate(., gene_symbol = factor(gene_symbol, levels = cdks_ccns)) %>%
  arrange(gene_symbol)
df_cdks_ccns_ht_tpm <- as.data.frame(df_cdks_ccns_ht_tpm)
rownames(df_cdks_ccns_ht_tpm) <- df_cdks_ccns_ht_tpm[["gene_symbol"]]
df_cdks_ccns_ht_tpm %<>%
  .[, -which(colnames(.) == "gene_symbol")]


# Setting the global parameter for heatmap title
ht_opt$TITLE_PADDING <- unit(c(8.5, 8.5), "points")

# Define output parameter
output_height <- 15
output_width <- 60

ht_tpm <- NULL
for (i in 1:length(cancers)) {
  target_cancer <- cancers[i]
  df_target_cancer_ht_tpm <- select(df_cdks_ccns_ht_tpm, starts_with(target_cancer))
  # cluster_cols: Calculate distance between each sample, and put the
  # result into df_target_cancer_ht_tpm_dist_column
  c_column <- dist(t(df_target_cancer_ht_tpm), method = "euclidean")
  df_target_cancer_ht_tpm_dist_column <- as.dist(c_column)
  # Execute clustering
  column_dend <- hclust(df_target_cancer_ht_tpm_dist_column, method = "average")
  # Define heatmap_param
  column_title <- paste0(target_cancer, "\n", "(", length(df_target_cancer_ht_tpm),
    ")")
  column_order <- colnames(df_target_cancer_ht_tpm)
  column_title_gp <- gpar(fontsize = 5, fontface = "bold", fill = i,
    col = "white", border = "white")
  row_title <- " "
  pre_ht_tpm <- Heatmap(as.matrix(df_target_cancer_ht_tpm), name = legend_title,
    col = col_fun, cluster_columns = column_dend, column_title = column_title,
    column_title_gp = column_title_gp, show_column_dend = FALSE, show_column_names = FALSE,
    row_split = 1:length(cdks_ccns), row_order = cdks_ccns, row_title = row_title,
    row_names_gp = gpar(fontsize = 12), row_names_side = "left", show_row_dend = FALSE,
    show_row_names = TRUE, height = unit(25, "cm"), heatmap_legend_param = heatmap_legend_param,
    use_raster = TRUE)
  ht_tpm <- ht_tpm + pre_ht_tpm
  pdf(paste0(top, "_CDKs", "_Clustered_Heatmap.pdf"), height = output_height,
    width = output_width)
  draw(ht_tpm)
  dev.off()
}

# Define cdks_ccns_list
cdks_ccns <- read_tsv("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/cdks.txt")[[1]]

# Reconstract dataframe
df_cdks_ccns_ht_tpm <- df_ht_tpm %>%
  filter_at(., vars(all_of("gene_symbol")), all_vars(. %in% cdks_ccns)) %>%
  mutate(., gene_symbol = factor(gene_symbol, levels = cdks_ccns)) %>%
  arrange(gene_symbol)
df_cdks_ccns_ht_tpm <- as.data.frame(df_cdks_ccns_ht_tpm)
rownames(df_cdks_ccns_ht_tpm) <- df_cdks_ccns_ht_tpm[["gene_symbol"]]
df_cdks_ccns_ht_tpm %<>%
  .[, -which(colnames(.) == "gene_symbol")]


# Setting the global parameter for heatmap title
ht_opt$TITLE_PADDING <- unit(c(8.5, 8.5), "points")

# Define output parameter
output_height <- 15
output_width <- 60

ht_tpm <- NULL
for (i in 1:length(cancers)) {
  target_cancer <- cancers[i]
  df_target_cancer_ht_tpm <- select(df_cdks_ccns_ht_tpm, starts_with(target_cancer))
  # cluster_cols: Calculate distance between each sample, and put the
  # result into df_target_cancer_ht_tpm_dist_column
  c_column <- dist(t(df_target_cancer_ht_tpm), method = "euclidean")
  df_target_cancer_ht_tpm_dist_column <- as.dist(c_column)
  # Execute clustering
  column_dend <- hclust(df_target_cancer_ht_tpm_dist_column, method = "average")
  # Define heatmap_param
  column_title <- paste0(target_cancer, "\n", "(", length(df_target_cancer_ht_tpm),
    ")")
  column_order <- colnames(df_target_cancer_ht_tpm)
  column_title_gp <- gpar(fontsize = 5, fontface = "bold", fill = i,
    col = "white", border = "white")
  row_title <- " "
  pre_ht_tpm <- Heatmap(as.matrix(df_target_cancer_ht_tpm), name = legend_title,
    col = col_fun, cluster_columns = column_dend, column_title = column_title,
    column_title_gp = column_title_gp, show_column_dend = FALSE, show_column_names = FALSE,
    row_split = 1:length(cdks_ccns), row_order = cdks_ccns, row_title = row_title,
    row_names_gp = gpar(fontsize = 11), row_names_side = "left", show_row_dend = FALSE,
    show_row_names = TRUE, height = unit(25, "cm"), heatmap_legend_param = heatmap_legend_param,
    use_raster = TRUE)
  ht_tpm <- ht_tpm + pre_ht_tpm
  pdf(paste0(top, "_CCNs", "_Clustered_Heatmap.pdf"), height = output_height,
    width = output_width)
  draw(ht_tpm)
  dev.off()
}

# Define legend_param
legend_title <- "gene expression (TPM)"
bottom <- 0
top <- 200
middle <- top/2
col_fun <- colorRamp2(c(bottom, middle, top), c("blue", "white", "red"))
heatmap_legend_param <- list(title_position = "leftcenter-rot", legend_height = unit(18,
  "cm"), grid_width = unit(1, "cm"), title_gp = gpar(fontsize = 20),
  labels_gp = gpar(fontsize = 20, font = 10))

# Define cdks_ccns_list
cdks_ccns <- read_tsv("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/cdks_ccns.txt")[[1]]

# Reconstract dataframe
df_cdks_ccns_ht_tpm <- df_ht_tpm %>%
  filter_at(., vars(all_of("gene_symbol")), all_vars(. %in% cdks_ccns)) %>%
  mutate(., gene_symbol = factor(gene_symbol, levels = cdks_ccns)) %>%
  arrange(gene_symbol)
df_cdks_ccns_ht_tpm <- as.data.frame(df_cdks_ccns_ht_tpm)
rownames(df_cdks_ccns_ht_tpm) <- df_cdks_ccns_ht_tpm[["gene_symbol"]]
df_cdks_ccns_ht_tpm %<>%
  .[, -which(colnames(.) == "gene_symbol")]

# Setting the global parameter for heatmap title
ht_opt$TITLE_PADDING <- unit(c(8.5, 8.5), "points")

# Define output parameter
output_height <- 25
output_width <- 45

ht_tpm <- NULL
for (i in 1:length(cancers)) {
  target_cancer <- cancers[i]
  df_target_cancer_ht_tpm <- select(df_cdks_ccns_ht_tpm, starts_with(target_cancer))
  # cluster_cols: Calculate distance between each sample, and put the
  # result into df_target_cancer_ht_tpm_dist_column
  c_column <- dist(t(df_target_cancer_ht_tpm), method = "euclidean")
  df_target_cancer_ht_tpm_dist_column <- as.dist(c_column)
  # Execute clustering
  column_dend <- hclust(df_target_cancer_ht_tpm_dist_column, method = "average")
  # Define heatmap_param
  column_title <- paste0(target_cancer, "\n", "(", length(df_target_cancer_ht_tpm),
    ")")
  column_order <- colnames(df_target_cancer_ht_tpm)
  column_title_gp <- gpar(fontsize = 5, fontface = "bold", fill = i,
    col = "white", border = "white")
  row_title <- " "
  pre_ht_tpm <- Heatmap(as.matrix(df_target_cancer_ht_tpm), name = legend_title,
    col = col_fun, cluster_columns = column_dend, column_title = column_title,
    column_title_gp = column_title_gp, show_column_dend = FALSE, show_column_names = FALSE,
    row_split = 1:length(cdks_ccns), row_order = cdks_ccns, row_title = row_title,
    row_names_gp = gpar(fontsize = 10), row_names_side = "left", show_row_dend = FALSE,
    show_row_names = TRUE, height = unit(50, "cm"), heatmap_legend_param = heatmap_legend_param,
    use_raster = TRUE)
  ht_tpm <- ht_tpm + pre_ht_tpm
  pdf(paste0(top, "_CDKs_CCNs", "_Clustered_Heatmap.pdf"), height = output_height,
    width = output_width)
  draw(ht_tpm)
  dev.off()
}

# Define cdks_ccns_list
cdks_ccns <- read_tsv("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/cdks.txt")[[1]]

# Reconstract dataframe
df_cdks_ccns_ht_tpm <- df_ht_tpm %>%
  filter_at(., vars(all_of("gene_symbol")), all_vars(. %in% cdks_ccns)) %>%
  mutate(., gene_symbol = factor(gene_symbol, levels = cdks_ccns)) %>%
  arrange(gene_symbol)
df_cdks_ccns_ht_tpm <- as.data.frame(df_cdks_ccns_ht_tpm)
rownames(df_cdks_ccns_ht_tpm) <- df_cdks_ccns_ht_tpm[["gene_symbol"]]
df_cdks_ccns_ht_tpm %<>%
  .[, -which(colnames(.) == "gene_symbol")]


# Setting the global parameter for heatmap title
ht_opt$TITLE_PADDING <- unit(c(8.5, 8.5), "points")

# Define output parameter
output_height <- 15
output_width <- 60

ht_tpm <- NULL
for (i in 1:length(cancers)) {
  target_cancer <- cancers[i]
  df_target_cancer_ht_tpm <- select(df_cdks_ccns_ht_tpm, starts_with(target_cancer))
  # cluster_cols: Calculate distance between each sample, and put the
  # result into df_target_cancer_ht_tpm_dist_column
  c_column <- dist(t(df_target_cancer_ht_tpm), method = "euclidean")
  df_target_cancer_ht_tpm_dist_column <- as.dist(c_column)
  # Execute clustering
  column_dend <- hclust(df_target_cancer_ht_tpm_dist_column, method = "average")
  # Define heatmap_param
  column_title <- paste0(target_cancer, "\n", "(", length(df_target_cancer_ht_tpm),
    ")")
  column_order <- colnames(df_target_cancer_ht_tpm)
  column_title_gp <- gpar(fontsize = 5, fontface = "bold", fill = i,
    col = "white", border = "white")
  row_title <- " "
  pre_ht_tpm <- Heatmap(as.matrix(df_target_cancer_ht_tpm), name = legend_title,
    col = col_fun, cluster_columns = column_dend, column_title = column_title,
    column_title_gp = column_title_gp, show_column_dend = FALSE, show_column_names = FALSE,
    row_split = 1:length(cdks_ccns), row_order = cdks_ccns, row_title = row_title,
    row_names_gp = gpar(fontsize = 12), row_names_side = "left", show_row_dend = FALSE,
    show_row_names = TRUE, height = unit(25, "cm"), heatmap_legend_param = heatmap_legend_param,
    use_raster = TRUE)
  ht_tpm <- ht_tpm + pre_ht_tpm
  pdf(paste0(top, "_CDKs", "_Clustered_Heatmap.pdf"), height = output_height,
    width = output_width)
  draw(ht_tpm)
  dev.off()
}

# Define cdks_ccns_list
cdks_ccns <- read_tsv("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/cdks.txt")[[1]]

# Reconstract dataframe
df_cdks_ccns_ht_tpm <- df_ht_tpm %>%
  filter_at(., vars(all_of("gene_symbol")), all_vars(. %in% cdks_ccns)) %>%
  mutate(., gene_symbol = factor(gene_symbol, levels = cdks_ccns)) %>%
  arrange(gene_symbol)
df_cdks_ccns_ht_tpm <- as.data.frame(df_cdks_ccns_ht_tpm)
rownames(df_cdks_ccns_ht_tpm) <- df_cdks_ccns_ht_tpm[["gene_symbol"]]
df_cdks_ccns_ht_tpm %<>%
  .[, -which(colnames(.) == "gene_symbol")]


# Setting the global parameter for heatmap title
ht_opt$TITLE_PADDING <- unit(c(8.5, 8.5), "points")

# Define output parameter
output_height <- 15
output_width <- 60

ht_tpm <- NULL
for (i in 1:length(cancers)) {
  target_cancer <- cancers[i]
  df_target_cancer_ht_tpm <- select(df_cdks_ccns_ht_tpm, starts_with(target_cancer))
  # cluster_cols: Calculate distance between each sample, and put the
  # result into df_target_cancer_ht_tpm_dist_column
  c_column <- dist(t(df_target_cancer_ht_tpm), method = "euclidean")
  df_target_cancer_ht_tpm_dist_column <- as.dist(c_column)
  # Execute clustering
  column_dend <- hclust(df_target_cancer_ht_tpm_dist_column, method = "average")
  # Define heatmap_param
  column_title <- paste0(target_cancer, "\n", "(", length(df_target_cancer_ht_tpm),
    ")")
  column_order <- colnames(df_target_cancer_ht_tpm)
  column_title_gp <- gpar(fontsize = 5, fontface = "bold", fill = i,
    col = "white", border = "white")
  row_title <- " "
  pre_ht_tpm <- Heatmap(as.matrix(df_target_cancer_ht_tpm), name = legend_title,
    col = col_fun, cluster_columns = column_dend, column_title = column_title,
    column_title_gp = column_title_gp, show_column_dend = FALSE, show_column_names = FALSE,
    row_split = 1:length(cdks_ccns), row_order = cdks_ccns, row_title = row_title,
    row_names_gp = gpar(fontsize = 11), row_names_side = "left", show_row_dend = FALSE,
    show_row_names = TRUE, height = unit(25, "cm"), heatmap_legend_param = heatmap_legend_param,
    use_raster = TRUE)
  ht_tpm <- ht_tpm + pre_ht_tpm
  pdf(paste0(top, "_CCNs", "_Clustered_Heatmap.pdf"), height = output_height,
    width = output_width)
  draw(ht_tpm)
  dev.off()
}

# Define legend_param
legend_title <- "gene expression (TPM)"
bottom <- 0
top <- 300
middle <- top/2
col_fun <- colorRamp2(c(bottom, middle, top), c("blue", "white", "red"))
heatmap_legend_param <- list(title_position = "leftcenter-rot", legend_height = unit(18,
  "cm"), grid_width = unit(1, "cm"), title_gp = gpar(fontsize = 20),
  labels_gp = gpar(fontsize = 20, font = 10))

# Define cdks_ccns_list
cdks_ccns <- read_tsv("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/cdks_ccns.txt")[[1]]

# Reconstract dataframe
df_cdks_ccns_ht_tpm <- df_ht_tpm %>%
  filter_at(., vars(all_of("gene_symbol")), all_vars(. %in% cdks_ccns)) %>%
  mutate(., gene_symbol = factor(gene_symbol, levels = cdks_ccns)) %>%
  arrange(gene_symbol)
df_cdks_ccns_ht_tpm <- as.data.frame(df_cdks_ccns_ht_tpm)
rownames(df_cdks_ccns_ht_tpm) <- df_cdks_ccns_ht_tpm[["gene_symbol"]]
df_cdks_ccns_ht_tpm %<>%
  .[, -which(colnames(.) == "gene_symbol")]

# Setting the global parameter for heatmap title
ht_opt$TITLE_PADDING <- unit(c(8.5, 8.5), "points")

# Define output parameter
output_height <- 25
output_width <- 45

ht_tpm <- NULL
for (i in 1:length(cancers)) {
  target_cancer <- cancers[i]
  df_target_cancer_ht_tpm <- select(df_cdks_ccns_ht_tpm, starts_with(target_cancer))
  # cluster_cols: Calculate distance between each sample, and put the
  # result into df_target_cancer_ht_tpm_dist_column
  c_column <- dist(t(df_target_cancer_ht_tpm), method = "euclidean")
  df_target_cancer_ht_tpm_dist_column <- as.dist(c_column)
  # Execute clustering
  column_dend <- hclust(df_target_cancer_ht_tpm_dist_column, method = "average")
  # Define heatmap_param
  column_title <- paste0(target_cancer, "\n", "(", length(df_target_cancer_ht_tpm),
    ")")
  column_order <- colnames(df_target_cancer_ht_tpm)
  column_title_gp <- gpar(fontsize = 5, fontface = "bold", fill = i,
    col = "white", border = "white")
  row_title <- " "
  pre_ht_tpm <- Heatmap(as.matrix(df_target_cancer_ht_tpm), name = legend_title,
    col = col_fun, cluster_columns = column_dend, column_title = column_title,
    column_title_gp = column_title_gp, show_column_dend = FALSE, show_column_names = FALSE,
    row_split = 1:length(cdks_ccns), row_order = cdks_ccns, row_title = row_title,
    row_names_gp = gpar(fontsize = 10), row_names_side = "left", show_row_dend = FALSE,
    show_row_names = TRUE, height = unit(50, "cm"), heatmap_legend_param = heatmap_legend_param,
    use_raster = TRUE)
  ht_tpm <- ht_tpm + pre_ht_tpm
  pdf(paste0(top, "_CDKs_CCNs", "_Clustered_Heatmap.pdf"), height = output_height,
    width = output_width)
  draw(ht_tpm)
  dev.off()
}

# Define cdks_ccns_list
cdks_ccns <- read_tsv("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/cdks.txt")[[1]]

# Reconstract dataframe
df_cdks_ccns_ht_tpm <- df_ht_tpm %>%
  filter_at(., vars(all_of("gene_symbol")), all_vars(. %in% cdks_ccns)) %>%
  mutate(., gene_symbol = factor(gene_symbol, levels = cdks_ccns)) %>%
  arrange(gene_symbol)
df_cdks_ccns_ht_tpm <- as.data.frame(df_cdks_ccns_ht_tpm)
rownames(df_cdks_ccns_ht_tpm) <- df_cdks_ccns_ht_tpm[["gene_symbol"]]
df_cdks_ccns_ht_tpm %<>%
  .[, -which(colnames(.) == "gene_symbol")]


# Setting the global parameter for heatmap title
ht_opt$TITLE_PADDING <- unit(c(8.5, 8.5), "points")

# Define output parameter
output_height <- 15
output_width <- 60

ht_tpm <- NULL
for (i in 1:length(cancers)) {
  target_cancer <- cancers[i]
  df_target_cancer_ht_tpm <- select(df_cdks_ccns_ht_tpm, starts_with(target_cancer))
  # cluster_cols: Calculate distance between each sample, and put the
  # result into df_target_cancer_ht_tpm_dist_column
  c_column <- dist(t(df_target_cancer_ht_tpm), method = "euclidean")
  df_target_cancer_ht_tpm_dist_column <- as.dist(c_column)
  # Execute clustering
  column_dend <- hclust(df_target_cancer_ht_tpm_dist_column, method = "average")
  # Define heatmap_param
  column_title <- paste0(target_cancer, "\n", "(", length(df_target_cancer_ht_tpm),
    ")")
  column_order <- colnames(df_target_cancer_ht_tpm)
  column_title_gp <- gpar(fontsize = 5, fontface = "bold", fill = i,
    col = "white", border = "white")
  row_title <- " "
  pre_ht_tpm <- Heatmap(as.matrix(df_target_cancer_ht_tpm), name = legend_title,
    col = col_fun, cluster_columns = column_dend, column_title = column_title,
    column_title_gp = column_title_gp, show_column_dend = FALSE, show_column_names = FALSE,
    row_split = 1:length(cdks_ccns), row_order = cdks_ccns, row_title = row_title,
    row_names_gp = gpar(fontsize = 12), row_names_side = "left", show_row_dend = FALSE,
    show_row_names = TRUE, height = unit(25, "cm"), heatmap_legend_param = heatmap_legend_param,
    use_raster = TRUE)
  ht_tpm <- ht_tpm + pre_ht_tpm
  pdf(paste0(top, "_CDKs", "_Clustered_Heatmap.pdf"), height = output_height,
    width = output_width)
  draw(ht_tpm)
  dev.off()
}

# Define cdks_ccns_list
cdks_ccns <- read_tsv("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/cdks.txt")[[1]]

# Reconstract dataframe
df_cdks_ccns_ht_tpm <- df_ht_tpm %>%
  filter_at(., vars(all_of("gene_symbol")), all_vars(. %in% cdks_ccns)) %>%
  mutate(., gene_symbol = factor(gene_symbol, levels = cdks_ccns)) %>%
  arrange(gene_symbol)
df_cdks_ccns_ht_tpm <- as.data.frame(df_cdks_ccns_ht_tpm)
rownames(df_cdks_ccns_ht_tpm) <- df_cdks_ccns_ht_tpm[["gene_symbol"]]
df_cdks_ccns_ht_tpm %<>%
  .[, -which(colnames(.) == "gene_symbol")]


# Setting the global parameter for heatmap title
ht_opt$TITLE_PADDING <- unit(c(8.5, 8.5), "points")

# Define output parameter
output_height <- 15
output_width <- 60

ht_tpm <- NULL
for (i in 1:length(cancers)) {
  target_cancer <- cancers[i]
  df_target_cancer_ht_tpm <- select(df_cdks_ccns_ht_tpm, starts_with(target_cancer))
  # cluster_cols: Calculate distance between each sample, and put the
  # result into df_target_cancer_ht_tpm_dist_column
  c_column <- dist(t(df_target_cancer_ht_tpm), method = "euclidean")
  df_target_cancer_ht_tpm_dist_column <- as.dist(c_column)
  # Execute clustering
  column_dend <- hclust(df_target_cancer_ht_tpm_dist_column, method = "average")
  # Define heatmap_param
  column_title <- paste0(target_cancer, "\n", "(", length(df_target_cancer_ht_tpm),
    ")")
  column_order <- colnames(df_target_cancer_ht_tpm)
  column_title_gp <- gpar(fontsize = 5, fontface = "bold", fill = i,
    col = "white", border = "white")
  row_title <- " "
  pre_ht_tpm <- Heatmap(as.matrix(df_target_cancer_ht_tpm), name = legend_title,
    col = col_fun, cluster_columns = column_dend, column_title = column_title,
    column_title_gp = column_title_gp, show_column_dend = FALSE, show_column_names = FALSE,
    row_split = 1:length(cdks_ccns), row_order = cdks_ccns, row_title = row_title,
    row_names_gp = gpar(fontsize = 11), row_names_side = "left", show_row_dend = FALSE,
    show_row_names = TRUE, height = unit(25, "cm"), heatmap_legend_param = heatmap_legend_param,
    use_raster = TRUE)
  ht_tpm <- ht_tpm + pre_ht_tpm
  pdf(paste0(top, "_CCNs", "_Clustered_Heatmap.pdf"), height = output_height,
    width = output_width)
  draw(ht_tpm)
  dev.off()
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Load package
library(tidyverse)
library(magrittr)
library(ComplexHeatmap)
library(circlize)

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

# Load data and select cdks_ccns samples
df_ht_tpm <- read_tsv("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/Datasets/lineup_cdks_ccns_tcga_tpm.txt")

# Define cancer_list
cancers <- unique(gsub("\\d", "", colnames(df_ht_tpm))) %>%
  .[-which(. == "gene_symbol")]

# Make Heatmap
# -------------------------------------------------------------------------
# Define Save directry
save_predir <- file.path("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer")
setwd(save_predir)
make_folder("Expression_Heatmap")
save_dir <- file.path(save_predir, "Expression_Heatmap")
setwd(save_dir)
make_folder("tpm")
setwd(file.path(save_dir, "tpm"))

# Ignores the warning message: Heatmap/annotation names are
# duplicated:Heatmap/annotation names are duplicated: gene expression
options(warn = -1)

# Define legend_param
legend_title <- "gene expression (TPM)"
bottom <- 0
top <- 50
middle <- top/2
col_fun <- colorRamp2(c(bottom, middle, top), c("blue", "white", "red"))
heatmap_legend_param <- list(title_position = "leftcenter-rot", legend_height = unit(18,
  "cm"), grid_width = unit(1, "cm"), title_gp = gpar(fontsize = 20),
  labels_gp = gpar(fontsize = 20, font = 10))

# Define cdks_ccns_list
cdks_ccns <- read_tsv("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/cdks_ccns.txt")[[1]]

# Reconstract dataframe
df_cdks_ccns_ht_tpm <- df_ht_tpm %>%
  filter_at(., vars(all_of("gene_symbol")), all_vars(. %in% cdks_ccns))
df_cdks_ccns_ht_tpm <- as.data.frame(df_cdks_ccns_ht_tpm)
rownames(df_cdks_ccns_ht_tpm) <- df_cdks_ccns_ht_tpm[["gene_symbol"]]
df_cdks_ccns_ht_tpm %<>%
  .[, -which(colnames(.) == "gene_symbol")]

# Setting the global parameter for heatmap title
ht_opt$TITLE_PADDING <- unit(c(8.5, 8.5), "points")

# Define output parameter
output_height <- 25
output_width <- 45

ht_tpm <- NULL
for (i in 1:length(cancers)) {
  target_cancer <- cancers[i]
  # Define heatmap_param
  df_target_cancer_ht_tpm <- select(df_cdks_ccns_ht_tpm, starts_with(target_cancer))
  column_title <- paste0(target_cancer, "\n", "(", length(df_target_cancer_ht_tpm),
    ")")
  column_order <- colnames(df_target_cancer_ht_tpm)
  column_title_gp <- gpar(fontsize = 12, fontface = "bold", fill = i,
    col = "white", border = "white")
  row_title <- " "
  pre_ht_tpm <- Heatmap(as.matrix(df_target_cancer_ht_tpm), name = legend_title,
    col = col_fun, column_order = column_order, column_title = column_title,
    column_title_gp = column_title_gp, show_column_dend = FALSE, show_column_names = FALSE,
    row_split = 1:length(cdks_ccns), row_order = cdks_ccns, row_title = row_title,
    row_names_gp = gpar(fontsize = 10), row_names_side = "left", show_row_dend = FALSE,
    show_row_names = TRUE, height = unit(50, "cm"), width = unit(3,
      "cm"), heatmap_legend_param = heatmap_legend_param, use_raster = TRUE)
  ht_tpm <- ht_tpm + pre_ht_tpm
  pdf(paste0(top, "_adj_CDKs_CCNs", "_ratio_Heatmap.pdf"), height = output_height,
    width = output_width)
  draw(ht_tpm)
  dev.off()
}

# Define cdks_ccns_list
cdks_ccns <- read_tsv("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/cdks.txt")[[1]]

# Reconstract dataframe
df_cdks_ccns_ht_tpm <- df_ht_tpm %>%
  filter_at(., vars(all_of("gene_symbol")), all_vars(. %in% cdks_ccns))
df_cdks_ccns_ht_tpm <- as.data.frame(df_cdks_ccns_ht_tpm)
rownames(df_cdks_ccns_ht_tpm) <- df_cdks_ccns_ht_tpm[["gene_symbol"]]
df_cdks_ccns_ht_tpm %<>%
  .[, -which(colnames(.) == "gene_symbol")]

# Setting the global parameter for heatmap title
ht_opt$TITLE_PADDING <- unit(c(8.5, 8.5), "points")

# Define output parameter
output_height <- 15
output_width <- 45

ht_tpm <- NULL
for (i in 1:length(cancers)) {
  target_cancer <- cancers[i]
  # Define heatmap_param
  df_target_cancer_ht_tpm <- select(df_cdks_ccns_ht_tpm, starts_with(target_cancer))
  column_title <- paste0(target_cancer, "\n", "(", length(df_target_cancer_ht_tpm),
    ")")
  column_order <- colnames(df_target_cancer_ht_tpm)
  column_title_gp <- gpar(fontsize = 12, fontface = "bold", fill = i,
    col = "white", border = "white")
  row_title <- " "
  pre_ht_tpm <- Heatmap(as.matrix(df_target_cancer_ht_tpm), name = legend_title,
    col = col_fun, column_order = column_order, column_title = column_title,
    column_title_gp = column_title_gp, show_column_dend = FALSE, show_column_names = FALSE,
    row_split = 1:length(cdks_ccns), row_order = cdks_ccns, row_title = row_title,
    row_names_gp = gpar(fontsize = 12), row_names_side = "left", show_row_dend = FALSE,
    show_row_names = TRUE, height = unit(25, "cm"), width = unit(3,
      "cm"), heatmap_legend_param = heatmap_legend_param, use_raster = TRUE)
  ht_tpm <- ht_tpm + pre_ht_tpm
  pdf(paste0(top, "_adj_CDKs", "_Heatmap.pdf"), height = output_height,
    width = output_width)
  draw(ht_tpm)
  dev.off()
}

# Define cdks_ccns_list
cdks_ccns <- read_tsv("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/ccns.txt")[[1]]

# Reconstract dataframe
df_cdks_ccns_ht_tpm <- df_ht_tpm %>%
  filter_at(., vars(all_of("gene_symbol")), all_vars(. %in% cdks_ccns))
df_cdks_ccns_ht_tpm <- as.data.frame(df_cdks_ccns_ht_tpm)
rownames(df_cdks_ccns_ht_tpm) <- df_cdks_ccns_ht_tpm[["gene_symbol"]]
df_cdks_ccns_ht_tpm %<>%
  .[, -which(colnames(.) == "gene_symbol")]

# Setting the global parameter for heatmap title
ht_opt$TITLE_PADDING <- unit(c(8.5, 8.5), "points")

# Define output parameter
output_height <- 15
output_width <- 45

ht_tpm <- NULL
for (i in 1:length(cancers)) {
  target_cancer <- cancers[i]
  # Define heatmap_param
  df_target_cancer_ht_tpm <- select(df_cdks_ccns_ht_tpm, starts_with(target_cancer))
  column_title <- paste0(target_cancer, "\n", "(", length(df_target_cancer_ht_tpm),
    ")")
  column_order <- colnames(df_target_cancer_ht_tpm)
  column_title_gp <- gpar(fontsize = 12, fontface = "bold", fill = i,
    col = "white", border = "white")
  row_title <- " "
  pre_ht_tpm <- Heatmap(as.matrix(df_target_cancer_ht_tpm), name = legend_title,
    col = col_fun, column_order = column_order, column_title = column_title,
    column_title_gp = column_title_gp, show_column_dend = FALSE, show_column_names = FALSE,
    row_split = 1:length(cdks_ccns), row_order = cdks_ccns, row_title = row_title,
    row_names_gp = gpar(fontsize = 11), row_names_side = "left", show_row_dend = FALSE,
    show_row_names = TRUE, height = unit(25, "cm"), width = unit(3,
      "cm"), heatmap_legend_param = heatmap_legend_param, use_raster = TRUE)
  ht_tpm <- ht_tpm + pre_ht_tpm
  pdf(paste0(top, "_adj_CCNs", "_Heatmap.pdf"), height = output_height,
    width = output_width)
  draw(ht_tpm)
  dev.off()
}

# Define legend_param
legend_title <- "gene expression (TPM)"
bottom <- 0
top <- 100
middle <- top/2
col_fun <- colorRamp2(c(bottom, middle, top), c("blue", "white", "red"))
heatmap_legend_param <- list(title_position = "leftcenter-rot", legend_height = unit(18,
  "cm"), grid_width = unit(1, "cm"), title_gp = gpar(fontsize = 20),
  labels_gp = gpar(fontsize = 20, font = 10))

# Define cdks_ccns_list
cdks_ccns <- read_tsv("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/cdks_ccns.txt")[[1]]

# Reconstract dataframe
df_cdks_ccns_ht_tpm <- df_ht_tpm %>%
  filter_at(., vars(all_of("gene_symbol")), all_vars(. %in% cdks_ccns))
df_cdks_ccns_ht_tpm <- as.data.frame(df_cdks_ccns_ht_tpm)
rownames(df_cdks_ccns_ht_tpm) <- df_cdks_ccns_ht_tpm[["gene_symbol"]]
df_cdks_ccns_ht_tpm %<>%
  .[, -which(colnames(.) == "gene_symbol")]

# Setting the global parameter for heatmap title
ht_opt$TITLE_PADDING <- unit(c(8.5, 8.5), "points")

# Define output parameter
output_height <- 25
output_width <- 45

ht_tpm <- NULL
for (i in 1:length(cancers)) {
  target_cancer <- cancers[i]
  # Define heatmap_param
  df_target_cancer_ht_tpm <- select(df_cdks_ccns_ht_tpm, starts_with(target_cancer))
  column_title <- paste0(target_cancer, "\n", "(", length(df_target_cancer_ht_tpm),
    ")")
  column_order <- colnames(df_target_cancer_ht_tpm)
  column_title_gp <- gpar(fontsize = 12, fontface = "bold", fill = i,
    col = "white", border = "white")
  row_title <- " "
  pre_ht_tpm <- Heatmap(as.matrix(df_target_cancer_ht_tpm), name = legend_title,
    col = col_fun, column_order = column_order, column_title = column_title,
    column_title_gp = column_title_gp, show_column_dend = FALSE, show_column_names = FALSE,
    row_split = 1:length(cdks_ccns), row_order = cdks_ccns, row_title = row_title,
    row_names_gp = gpar(fontsize = 10), row_names_side = "left", show_row_dend = FALSE,
    show_row_names = TRUE, height = unit(50, "cm"), width = unit(3,
      "cm"), heatmap_legend_param = heatmap_legend_param, use_raster = TRUE)
  ht_tpm <- ht_tpm + pre_ht_tpm
  pdf(paste0(top, "_adj_CDKs_CCNs", "_ratio_Heatmap.pdf"), height = output_height,
    width = output_width)
  draw(ht_tpm)
  dev.off()
}

# Define cdks_ccns_list
cdks_ccns <- read_tsv("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/cdks.txt")[[1]]

# Reconstract dataframe
df_cdks_ccns_ht_tpm <- df_ht_tpm %>%
  filter_at(., vars(all_of("gene_symbol")), all_vars(. %in% cdks_ccns))
df_cdks_ccns_ht_tpm <- as.data.frame(df_cdks_ccns_ht_tpm)
rownames(df_cdks_ccns_ht_tpm) <- df_cdks_ccns_ht_tpm[["gene_symbol"]]
df_cdks_ccns_ht_tpm %<>%
  .[, -which(colnames(.) == "gene_symbol")]

# Setting the global parameter for heatmap title
ht_opt$TITLE_PADDING <- unit(c(8.5, 8.5), "points")

# Define output parameter
output_height <- 15
output_width <- 45

ht_tpm <- NULL
for (i in 1:length(cancers)) {
  target_cancer <- cancers[i]
  # Define heatmap_param
  df_target_cancer_ht_tpm <- select(df_cdks_ccns_ht_tpm, starts_with(target_cancer))
  column_title <- paste0(target_cancer, "\n", "(", length(df_target_cancer_ht_tpm),
    ")")
  column_order <- colnames(df_target_cancer_ht_tpm)
  column_title_gp <- gpar(fontsize = 12, fontface = "bold", fill = i,
    col = "white", border = "white")
  row_title <- " "
  pre_ht_tpm <- Heatmap(as.matrix(df_target_cancer_ht_tpm), name = legend_title,
    col = col_fun, column_order = column_order, column_title = column_title,
    column_title_gp = column_title_gp, show_column_dend = FALSE, show_column_names = FALSE,
    row_split = 1:length(cdks_ccns), row_order = cdks_ccns, row_title = row_title,
    row_names_gp = gpar(fontsize = 12), row_names_side = "left", show_row_dend = FALSE,
    show_row_names = TRUE, height = unit(25, "cm"), width = unit(3,
      "cm"), heatmap_legend_param = heatmap_legend_param, use_raster = TRUE)
  ht_tpm <- ht_tpm + pre_ht_tpm
  pdf(paste0(top, "_adj_CDKs", "_Heatmap.pdf"), height = output_height,
    width = output_width)
  draw(ht_tpm)
  dev.off()
}

# Define cdks_ccns_list
cdks_ccns <- read_tsv("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/ccns.txt")[[1]]

# Reconstract dataframe
df_cdks_ccns_ht_tpm <- df_ht_tpm %>%
  filter_at(., vars(all_of("gene_symbol")), all_vars(. %in% cdks_ccns))
df_cdks_ccns_ht_tpm <- as.data.frame(df_cdks_ccns_ht_tpm)
rownames(df_cdks_ccns_ht_tpm) <- df_cdks_ccns_ht_tpm[["gene_symbol"]]
df_cdks_ccns_ht_tpm %<>%
  .[, -which(colnames(.) == "gene_symbol")]

# Setting the global parameter for heatmap title
ht_opt$TITLE_PADDING <- unit(c(8.5, 8.5), "points")

# Define output parameter
output_height <- 15
output_width <- 45

ht_tpm <- NULL
for (i in 1:length(cancers)) {
  target_cancer <- cancers[i]
  # Define heatmap_param
  df_target_cancer_ht_tpm <- select(df_cdks_ccns_ht_tpm, starts_with(target_cancer))
  column_title <- paste0(target_cancer, "\n", "(", length(df_target_cancer_ht_tpm),
    ")")
  column_order <- colnames(df_target_cancer_ht_tpm)
  column_title_gp <- gpar(fontsize = 12, fontface = "bold", fill = i,
    col = "white", border = "white")
  row_title <- " "
  pre_ht_tpm <- Heatmap(as.matrix(df_target_cancer_ht_tpm), name = legend_title,
    col = col_fun, column_order = column_order, column_title = column_title,
    column_title_gp = column_title_gp, show_column_dend = FALSE, show_column_names = FALSE,
    row_split = 1:length(cdks_ccns), row_order = cdks_ccns, row_title = row_title,
    row_names_gp = gpar(fontsize = 11), row_names_side = "left", show_row_dend = FALSE,
    show_row_names = TRUE, height = unit(25, "cm"), width = unit(3,
      "cm"), heatmap_legend_param = heatmap_legend_param, use_raster = TRUE)
  ht_tpm <- ht_tpm + pre_ht_tpm
  pdf(paste0(top, "_adj_CCNs", "_Heatmap.pdf"), height = output_height,
    width = output_width)
  draw(ht_tpm)
  dev.off()
}

# Define legend_param
legend_title <- "gene expression (TPM)"
bottom <- 0
top <- 150
middle <- top/2
col_fun <- colorRamp2(c(bottom, middle, top), c("blue", "white", "red"))
heatmap_legend_param <- list(title_position = "leftcenter-rot", legend_height = unit(18,
  "cm"), grid_width = unit(1, "cm"), title_gp = gpar(fontsize = 20),
  labels_gp = gpar(fontsize = 20, font = 10))

# Define cdks_ccns_list
cdks_ccns <- read_tsv("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/cdks_ccns.txt")[[1]]

# Reconstract dataframe
df_cdks_ccns_ht_tpm <- df_ht_tpm %>%
  filter_at(., vars(all_of("gene_symbol")), all_vars(. %in% cdks_ccns))
df_cdks_ccns_ht_tpm <- as.data.frame(df_cdks_ccns_ht_tpm)
rownames(df_cdks_ccns_ht_tpm) <- df_cdks_ccns_ht_tpm[["gene_symbol"]]
df_cdks_ccns_ht_tpm %<>%
  .[, -which(colnames(.) == "gene_symbol")]

# Setting the global parameter for heatmap title
ht_opt$TITLE_PADDING <- unit(c(8.5, 8.5), "points")

# Define output parameter
output_height <- 25
output_width <- 45

ht_tpm <- NULL
for (i in 1:length(cancers)) {
  target_cancer <- cancers[i]
  # Define heatmap_param
  df_target_cancer_ht_tpm <- select(df_cdks_ccns_ht_tpm, starts_with(target_cancer))
  column_title <- paste0(target_cancer, "\n", "(", length(df_target_cancer_ht_tpm),
    ")")
  column_order <- colnames(df_target_cancer_ht_tpm)
  column_title_gp <- gpar(fontsize = 12, fontface = "bold", fill = i,
    col = "white", border = "white")
  row_title <- " "
  pre_ht_tpm <- Heatmap(as.matrix(df_target_cancer_ht_tpm), name = legend_title,
    col = col_fun, column_order = column_order, column_title = column_title,
    column_title_gp = column_title_gp, show_column_dend = FALSE, show_column_names = FALSE,
    row_split = 1:length(cdks_ccns), row_order = cdks_ccns, row_title = row_title,
    row_names_gp = gpar(fontsize = 10), row_names_side = "left", show_row_dend = FALSE,
    show_row_names = TRUE, height = unit(50, "cm"), width = unit(3,
      "cm"), heatmap_legend_param = heatmap_legend_param, use_raster = TRUE)
  ht_tpm <- ht_tpm + pre_ht_tpm
  pdf(paste0(top, "_adj_CDKs_CCNs", "_ratio_Heatmap.pdf"), height = output_height,
    width = output_width)
  draw(ht_tpm)
  dev.off()
}

# Define cdks_ccns_list
cdks_ccns <- read_tsv("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/cdks.txt")[[1]]

# Reconstract dataframe
df_cdks_ccns_ht_tpm <- df_ht_tpm %>%
  filter_at(., vars(all_of("gene_symbol")), all_vars(. %in% cdks_ccns))
df_cdks_ccns_ht_tpm <- as.data.frame(df_cdks_ccns_ht_tpm)
rownames(df_cdks_ccns_ht_tpm) <- df_cdks_ccns_ht_tpm[["gene_symbol"]]
df_cdks_ccns_ht_tpm %<>%
  .[, -which(colnames(.) == "gene_symbol")]

# Setting the global parameter for heatmap title
ht_opt$TITLE_PADDING <- unit(c(8.5, 8.5), "points")

# Define output parameter
output_height <- 15
output_width <- 45

ht_tpm <- NULL
for (i in 1:length(cancers)) {
  target_cancer <- cancers[i]
  # Define heatmap_param
  df_target_cancer_ht_tpm <- select(df_cdks_ccns_ht_tpm, starts_with(target_cancer))
  column_title <- paste0(target_cancer, "\n", "(", length(df_target_cancer_ht_tpm),
    ")")
  column_order <- colnames(df_target_cancer_ht_tpm)
  column_title_gp <- gpar(fontsize = 12, fontface = "bold", fill = i,
    col = "white", border = "white")
  row_title <- " "
  pre_ht_tpm <- Heatmap(as.matrix(df_target_cancer_ht_tpm), name = legend_title,
    col = col_fun, column_order = column_order, column_title = column_title,
    column_title_gp = column_title_gp, show_column_dend = FALSE, show_column_names = FALSE,
    row_split = 1:length(cdks_ccns), row_order = cdks_ccns, row_title = row_title,
    row_names_gp = gpar(fontsize = 12), row_names_side = "left", show_row_dend = FALSE,
    show_row_names = TRUE, height = unit(25, "cm"), width = unit(3,
      "cm"), heatmap_legend_param = heatmap_legend_param, use_raster = TRUE)
  ht_tpm <- ht_tpm + pre_ht_tpm
  pdf(paste0(top, "_adj_CDKs", "_Heatmap.pdf"), height = output_height,
    width = output_width)
  draw(ht_tpm)
  dev.off()
}

# Define cdks_ccns_list
cdks_ccns <- read_tsv("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/ccns.txt")[[1]]

# Reconstract dataframe
df_cdks_ccns_ht_tpm <- df_ht_tpm %>%
  filter_at(., vars(all_of("gene_symbol")), all_vars(. %in% cdks_ccns))
df_cdks_ccns_ht_tpm <- as.data.frame(df_cdks_ccns_ht_tpm)
rownames(df_cdks_ccns_ht_tpm) <- df_cdks_ccns_ht_tpm[["gene_symbol"]]
df_cdks_ccns_ht_tpm %<>%
  .[, -which(colnames(.) == "gene_symbol")]

# Setting the global parameter for heatmap title
ht_opt$TITLE_PADDING <- unit(c(8.5, 8.5), "points")

# Define output parameter
output_height <- 15
output_width <- 45

ht_tpm <- NULL
for (i in 1:length(cancers)) {
  target_cancer <- cancers[i]
  # Define heatmap_param
  df_target_cancer_ht_tpm <- select(df_cdks_ccns_ht_tpm, starts_with(target_cancer))
  column_title <- paste0(target_cancer, "\n", "(", length(df_target_cancer_ht_tpm),
    ")")
  column_order <- colnames(df_target_cancer_ht_tpm)
  column_title_gp <- gpar(fontsize = 12, fontface = "bold", fill = i,
    col = "white", border = "white")
  row_title <- " "
  pre_ht_tpm <- Heatmap(as.matrix(df_target_cancer_ht_tpm), name = legend_title,
    col = col_fun, column_order = column_order, column_title = column_title,
    column_title_gp = column_title_gp, show_column_dend = FALSE, show_column_names = FALSE,
    row_split = 1:length(cdks_ccns), row_order = cdks_ccns, row_title = row_title,
    row_names_gp = gpar(fontsize = 11), row_names_side = "left", show_row_dend = FALSE,
    show_row_names = TRUE, height = unit(25, "cm"), width = unit(3,
      "cm"), heatmap_legend_param = heatmap_legend_param, use_raster = TRUE)
  ht_tpm <- ht_tpm + pre_ht_tpm
  pdf(paste0(top, "_adj_CCNs", "_Heatmap.pdf"), height = output_height,
    width = output_width)
  draw(ht_tpm)
  dev.off()
}

# Define legend_param
legend_title <- "gene expression (TPM)"
bottom <- 0
top <- 200
middle <- top/2
col_fun <- colorRamp2(c(bottom, middle, top), c("blue", "white", "red"))
heatmap_legend_param <- list(title_position = "leftcenter-rot", legend_height = unit(18,
  "cm"), grid_width = unit(1, "cm"), title_gp = gpar(fontsize = 20),
  labels_gp = gpar(fontsize = 20, font = 10))

# Define cdks_ccns_list
cdks_ccns <- read_tsv("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/cdks_ccns.txt")[[1]]

# Reconstract dataframe
df_cdks_ccns_ht_tpm <- df_ht_tpm %>%
  filter_at(., vars(all_of("gene_symbol")), all_vars(. %in% cdks_ccns))
df_cdks_ccns_ht_tpm <- as.data.frame(df_cdks_ccns_ht_tpm)
rownames(df_cdks_ccns_ht_tpm) <- df_cdks_ccns_ht_tpm[["gene_symbol"]]
df_cdks_ccns_ht_tpm %<>%
  .[, -which(colnames(.) == "gene_symbol")]

# Setting the global parameter for heatmap title
ht_opt$TITLE_PADDING <- unit(c(8.5, 8.5), "points")

# Define output parameter
output_height <- 25
output_width <- 45

ht_tpm <- NULL
for (i in 1:length(cancers)) {
  target_cancer <- cancers[i]
  # Define heatmap_param
  df_target_cancer_ht_tpm <- select(df_cdks_ccns_ht_tpm, starts_with(target_cancer))
  column_title <- paste0(target_cancer, "\n", "(", length(df_target_cancer_ht_tpm),
    ")")
  column_order <- colnames(df_target_cancer_ht_tpm)
  column_title_gp <- gpar(fontsize = 12, fontface = "bold", fill = i,
    col = "white", border = "white")
  row_title <- " "
  pre_ht_tpm <- Heatmap(as.matrix(df_target_cancer_ht_tpm), name = legend_title,
    col = col_fun, column_order = column_order, column_title = column_title,
    column_title_gp = column_title_gp, show_column_dend = FALSE, show_column_names = FALSE,
    row_split = 1:length(cdks_ccns), row_order = cdks_ccns, row_title = row_title,
    row_names_gp = gpar(fontsize = 10), row_names_side = "left", show_row_dend = FALSE,
    show_row_names = TRUE, height = unit(50, "cm"), width = unit(3,
      "cm"), heatmap_legend_param = heatmap_legend_param, use_raster = TRUE)
  ht_tpm <- ht_tpm + pre_ht_tpm
  pdf(paste0(top, "_adj_CDKs_CCNs", "_ratio_Heatmap.pdf"), height = output_height,
    width = output_width)
  draw(ht_tpm)
  dev.off()
}

# Define cdks_ccns_list
cdks_ccns <- read_tsv("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/cdks.txt")[[1]]

# Reconstract dataframe
df_cdks_ccns_ht_tpm <- df_ht_tpm %>%
  filter_at(., vars(all_of("gene_symbol")), all_vars(. %in% cdks_ccns))
df_cdks_ccns_ht_tpm <- as.data.frame(df_cdks_ccns_ht_tpm)
rownames(df_cdks_ccns_ht_tpm) <- df_cdks_ccns_ht_tpm[["gene_symbol"]]
df_cdks_ccns_ht_tpm %<>%
  .[, -which(colnames(.) == "gene_symbol")]

# Setting the global parameter for heatmap title
ht_opt$TITLE_PADDING <- unit(c(8.5, 8.5), "points")

# Define output parameter
output_height <- 15
output_width <- 45

ht_tpm <- NULL
for (i in 1:length(cancers)) {
  target_cancer <- cancers[i]
  # Define heatmap_param
  df_target_cancer_ht_tpm <- select(df_cdks_ccns_ht_tpm, starts_with(target_cancer))
  column_title <- paste0(target_cancer, "\n", "(", length(df_target_cancer_ht_tpm),
    ")")
  column_order <- colnames(df_target_cancer_ht_tpm)
  column_title_gp <- gpar(fontsize = 12, fontface = "bold", fill = i,
    col = "white", border = "white")
  row_title <- " "
  pre_ht_tpm <- Heatmap(as.matrix(df_target_cancer_ht_tpm), name = legend_title,
    col = col_fun, column_order = column_order, column_title = column_title,
    column_title_gp = column_title_gp, show_column_dend = FALSE, show_column_names = FALSE,
    row_split = 1:length(cdks_ccns), row_order = cdks_ccns, row_title = row_title,
    row_names_gp = gpar(fontsize = 12), row_names_side = "left", show_row_dend = FALSE,
    show_row_names = TRUE, height = unit(25, "cm"), width = unit(3,
      "cm"), heatmap_legend_param = heatmap_legend_param, use_raster = TRUE)
  ht_tpm <- ht_tpm + pre_ht_tpm
  pdf(paste0(top, "_adj_CDKs", "_Heatmap.pdf"), height = output_height,
    width = output_width)
  draw(ht_tpm)
  dev.off()
}

# Define cdks_ccns_list
cdks_ccns <- read_tsv("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/ccns.txt")[[1]]

# Reconstract dataframe
df_cdks_ccns_ht_tpm <- df_ht_tpm %>%
  filter_at(., vars(all_of("gene_symbol")), all_vars(. %in% cdks_ccns))
df_cdks_ccns_ht_tpm <- as.data.frame(df_cdks_ccns_ht_tpm)
rownames(df_cdks_ccns_ht_tpm) <- df_cdks_ccns_ht_tpm[["gene_symbol"]]
df_cdks_ccns_ht_tpm %<>%
  .[, -which(colnames(.) == "gene_symbol")]

# Setting the global parameter for heatmap title
ht_opt$TITLE_PADDING <- unit(c(8.5, 8.5), "points")

# Define output parameter
output_height <- 15
output_width <- 45

ht_tpm <- NULL
for (i in 1:length(cancers)) {
  target_cancer <- cancers[i]
  # Define heatmap_param
  df_target_cancer_ht_tpm <- select(df_cdks_ccns_ht_tpm, starts_with(target_cancer))
  column_title <- paste0(target_cancer, "\n", "(", length(df_target_cancer_ht_tpm),
    ")")
  column_order <- colnames(df_target_cancer_ht_tpm)
  column_title_gp <- gpar(fontsize = 12, fontface = "bold", fill = i,
    col = "white", border = "white")
  row_title <- " "
  pre_ht_tpm <- Heatmap(as.matrix(df_target_cancer_ht_tpm), name = legend_title,
    col = col_fun, column_order = column_order, column_title = column_title,
    column_title_gp = column_title_gp, show_column_dend = FALSE, show_column_names = FALSE,
    row_split = 1:length(cdks_ccns), row_order = cdks_ccns, row_title = row_title,
    row_names_gp = gpar(fontsize = 11), row_names_side = "left", show_row_dend = FALSE,
    show_row_names = TRUE, height = unit(25, "cm"), width = unit(3,
      "cm"), heatmap_legend_param = heatmap_legend_param, use_raster = TRUE)
  ht_tpm <- ht_tpm + pre_ht_tpm
  pdf(paste0(top, "_adj_CCNs", "_Heatmap.pdf"), height = output_height,
    width = output_width)
  draw(ht_tpm)
  dev.off()
}

# Define legend_param
legend_title <- "gene expression (TPM)"
bottom <- 0
top <- 300
middle <- top/2
col_fun <- colorRamp2(c(bottom, middle, top), c("blue", "white", "red"))
heatmap_legend_param <- list(title_position = "leftcenter-rot", legend_height = unit(18,
  "cm"), grid_width = unit(1, "cm"), title_gp = gpar(fontsize = 20),
  labels_gp = gpar(fontsize = 20, font = 10))

# Define cdks_ccns_list
cdks_ccns <- read_tsv("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/cdks_ccns.txt")[[1]]

# Reconstract dataframe
df_cdks_ccns_ht_tpm <- df_ht_tpm %>%
  filter_at(., vars(all_of("gene_symbol")), all_vars(. %in% cdks_ccns))
df_cdks_ccns_ht_tpm <- as.data.frame(df_cdks_ccns_ht_tpm)
rownames(df_cdks_ccns_ht_tpm) <- df_cdks_ccns_ht_tpm[["gene_symbol"]]
df_cdks_ccns_ht_tpm %<>%
  .[, -which(colnames(.) == "gene_symbol")]

# Setting the global parameter for heatmap title
ht_opt$TITLE_PADDING <- unit(c(8.5, 8.5), "points")

# Define output parameter
output_height <- 25
output_width <- 45

ht_tpm <- NULL
for (i in 1:length(cancers)) {
  target_cancer <- cancers[i]
  # Define heatmap_param
  df_target_cancer_ht_tpm <- select(df_cdks_ccns_ht_tpm, starts_with(target_cancer))
  column_title <- paste0(target_cancer, "\n", "(", length(df_target_cancer_ht_tpm),
    ")")
  column_order <- colnames(df_target_cancer_ht_tpm)
  column_title_gp <- gpar(fontsize = 12, fontface = "bold", fill = i,
    col = "white", border = "white")
  row_title <- " "
  pre_ht_tpm <- Heatmap(as.matrix(df_target_cancer_ht_tpm), name = legend_title,
    col = col_fun, column_order = column_order, column_title = column_title,
    column_title_gp = column_title_gp, show_column_dend = FALSE, show_column_names = FALSE,
    row_split = 1:length(cdks_ccns), row_order = cdks_ccns, row_title = row_title,
    row_names_gp = gpar(fontsize = 10), row_names_side = "left", show_row_dend = FALSE,
    show_row_names = TRUE, height = unit(50, "cm"), width = unit(3,
      "cm"), heatmap_legend_param = heatmap_legend_param, use_raster = TRUE)
  ht_tpm <- ht_tpm + pre_ht_tpm
  pdf(paste0(top, "_adj_CDKs_CCNs", "_ratio_Heatmap.pdf"), height = output_height,
    width = output_width)
  draw(ht_tpm)
  dev.off()
}

# Define cdks_ccns_list
cdks_ccns <- read_tsv("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/cdks.txt")[[1]]

# Reconstract dataframe
df_cdks_ccns_ht_tpm <- df_ht_tpm %>%
  filter_at(., vars(all_of("gene_symbol")), all_vars(. %in% cdks_ccns))
df_cdks_ccns_ht_tpm <- as.data.frame(df_cdks_ccns_ht_tpm)
rownames(df_cdks_ccns_ht_tpm) <- df_cdks_ccns_ht_tpm[["gene_symbol"]]
df_cdks_ccns_ht_tpm %<>%
  .[, -which(colnames(.) == "gene_symbol")]

# Setting the global parameter for heatmap title
ht_opt$TITLE_PADDING <- unit(c(8.5, 8.5), "points")

# Define output parameter
output_height <- 15
output_width <- 45

ht_tpm <- NULL
for (i in 1:length(cancers)) {
  target_cancer <- cancers[i]
  # Define heatmap_param
  df_target_cancer_ht_tpm <- select(df_cdks_ccns_ht_tpm, starts_with(target_cancer))
  column_title <- paste0(target_cancer, "\n", "(", length(df_target_cancer_ht_tpm),
    ")")
  column_order <- colnames(df_target_cancer_ht_tpm)
  column_title_gp <- gpar(fontsize = 12, fontface = "bold", fill = i,
    col = "white", border = "white")
  row_title <- " "
  pre_ht_tpm <- Heatmap(as.matrix(df_target_cancer_ht_tpm), name = legend_title,
    col = col_fun, column_order = column_order, column_title = column_title,
    column_title_gp = column_title_gp, show_column_dend = FALSE, show_column_names = FALSE,
    row_split = 1:length(cdks_ccns), row_order = cdks_ccns, row_title = row_title,
    row_names_gp = gpar(fontsize = 12), row_names_side = "left", show_row_dend = FALSE,
    show_row_names = TRUE, height = unit(25, "cm"), width = unit(3,
      "cm"), heatmap_legend_param = heatmap_legend_param, use_raster = TRUE)
  ht_tpm <- ht_tpm + pre_ht_tpm
  pdf(paste0(top, "_adj_CDKs", "_Heatmap.pdf"), height = output_height,
    width = output_width)
  draw(ht_tpm)
  dev.off()
}

# Define cdks_ccns_list
cdks_ccns <- read_tsv("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/ccns.txt")[[1]]

# Reconstract dataframe
df_cdks_ccns_ht_tpm <- df_ht_tpm %>%
  filter_at(., vars(all_of("gene_symbol")), all_vars(. %in% cdks_ccns))
df_cdks_ccns_ht_tpm <- as.data.frame(df_cdks_ccns_ht_tpm)
rownames(df_cdks_ccns_ht_tpm) <- df_cdks_ccns_ht_tpm[["gene_symbol"]]
df_cdks_ccns_ht_tpm %<>%
  .[, -which(colnames(.) == "gene_symbol")]

# Setting the global parameter for heatmap title
ht_opt$TITLE_PADDING <- unit(c(8.5, 8.5), "points")

# Define output parameter
output_height <- 15
output_width <- 45

ht_tpm <- NULL
for (i in 1:length(cancers)) {
  target_cancer <- cancers[i]
  # Define heatmap_param
  df_target_cancer_ht_tpm <- select(df_cdks_ccns_ht_tpm, starts_with(target_cancer))
  column_title <- paste0(target_cancer, "\n", "(", length(df_target_cancer_ht_tpm),
    ")")
  column_order <- colnames(df_target_cancer_ht_tpm)
  column_title_gp <- gpar(fontsize = 12, fontface = "bold", fill = i,
    col = "white", border = "white")
  row_title <- " "
  pre_ht_tpm <- Heatmap(as.matrix(df_target_cancer_ht_tpm), name = legend_title,
    col = col_fun, column_order = column_order, column_title = column_title,
    column_title_gp = column_title_gp, show_column_dend = FALSE, show_column_names = FALSE,
    row_split = 1:length(cdks_ccns), row_order = cdks_ccns, row_title = row_title,
    row_names_gp = gpar(fontsize = 11), row_names_side = "left", show_row_dend = FALSE,
    show_row_names = TRUE, height = unit(25, "cm"), width = unit(3,
      "cm"), heatmap_legend_param = heatmap_legend_param, use_raster = TRUE)
  ht_tpm <- ht_tpm + pre_ht_tpm
  pdf(paste0(top, "_adj_CCNs", "_Heatmap.pdf"), height = output_height,
    width = output_width)
  draw(ht_tpm)
  dev.off()
}

traceback()



