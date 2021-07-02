## We need to revise data sets from recount3 because TCGA datasets are including 
## a lot of variable samples. Those data are including FFPE samples, duplicates, etc.
## Before analysis, We need to curate samples. Fortunately, TCGA group published some articles using TCGA datasets,
## I curated samples based on their data.

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Load package
library(tidyverse)
library(magrittr)
library(readxl)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "acc"
df_reference <- "table_s1.xlsx"
sheet_reference <- "Data Overview"
skip_row <- 6
coln_reference <- "TCGA_ID"
except_row <- "Platform"
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "tcga.gdc_cases.submitter_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Remove all duplicate samples because we cannot confirm the reason
# why they selected one of the duplicate samples
lst_tcga_barcode <- meta_recount3[["tcga.gdc_cases.submitter_id"]]
dupl_barcode <- lst_tcga_barcode[duplicated(lst_tcga_barcode)]
if (length(dupl_barcode) != 0) {
  nondupl_barcode <- lst_tcga_barcode[-which(lst_tcga_barcode %in% dupl_barcode)]
} else {
  nondupl_barcode <- lst_tcga_barcode
}
meta_recount3 %<>%
  filter_at(., vars(all_of("tcga.gdc_cases.submitter_id")), all_vars(. %in% nondupl_barcode))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Change the colnames from external_id to tcga_barcode
lst_external_id <- colnames(counts_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  counts_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

lst_tpm_external_id <- colnames(tpm_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  tpm_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "blca_mibc"
df_reference <- "table_s1.xlsx"
sheet_reference <- "Master table"
skip_row <- 0
coln_reference <- "BCR patient uuid"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- tolower(meta_reference[[coln_reference]])

# Define recount3 parameters and load meta_data
target_recount3 <- "blca"
coln_recount3 <- "tcga.gdc_cases.case_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Remove all duplicate samples because we cannot confirm the reason
# why they selected one of the duplicate samples
lst_tcga_barcode <- meta_recount3[["tcga.gdc_cases.submitter_id"]]
dupl_barcode <- lst_tcga_barcode[duplicated(lst_tcga_barcode)]
if (length(dupl_barcode) != 0) {
  nondupl_barcode <- lst_tcga_barcode[-which(lst_tcga_barcode %in% dupl_barcode)]
} else {
  nondupl_barcode <- lst_tcga_barcode
}
meta_recount3 %<>%
  filter_at(., vars(all_of("tcga.gdc_cases.submitter_id")), all_vars(. %in% nondupl_barcode))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Change the colnames from external_id to tcga_barcode
lst_external_id <- colnames(counts_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  counts_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

lst_tpm_external_id <- colnames(tpm_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  tpm_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "brca"
df_reference <- "table_s1.xlsx"
sheet_reference <- "Suppl. Table 1"
skip_row <- 2
coln_reference <- "mRNA"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "tcga.tcga_barcode"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Remove all duplicate samples because we cannot confirm the reason
# why they selected one of the duplicate samples
lst_tcga_barcode <- meta_recount3[["tcga.gdc_cases.submitter_id"]]
dupl_barcode <- lst_tcga_barcode[duplicated(lst_tcga_barcode)]
if (length(dupl_barcode) != 0) {
  nondupl_barcode <- lst_tcga_barcode[-which(lst_tcga_barcode %in% dupl_barcode)]
} else {
  nondupl_barcode <- lst_tcga_barcode
}
meta_recount3 %<>%
  filter_at(., vars(all_of("tcga.gdc_cases.submitter_id")), all_vars(. %in% nondupl_barcode))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Change the colnames from external_id to tcga_barcode
lst_external_id <- colnames(counts_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  counts_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

lst_tpm_external_id <- colnames(tpm_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  tpm_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "chol"
df_reference <- "tables.xlsx"
sheet_reference <- "ST1. CHOL Sample Info"
skip_row <- 3
coln_reference <- "...1"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
meta_reference %<>% filter_at(., vars(all_of("Sample Type")), all_vars(. == "Primary solid Tumor"))
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "tcga.gdc_cases.submitter_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Remove all duplicate samples because we cannot confirm the reason
# why they selected one of the duplicate samples
lst_tcga_barcode <- meta_recount3[["tcga.gdc_cases.submitter_id"]]
dupl_barcode <- lst_tcga_barcode[duplicated(lst_tcga_barcode)]
if (length(dupl_barcode) != 0) {
  nondupl_barcode <- lst_tcga_barcode[-which(lst_tcga_barcode %in% dupl_barcode)]
} else {
  nondupl_barcode <- lst_tcga_barcode
}
meta_recount3 %<>%
  filter_at(., vars(all_of("tcga.gdc_cases.submitter_id")), all_vars(. %in% nondupl_barcode))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Change the colnames from external_id to tcga_barcode
lst_external_id <- colnames(counts_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  counts_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

lst_tpm_external_id <- colnames(tpm_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  tpm_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    "nature11252-s2", dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "coad_read"
df_reference <- "2011-11-14592C-Sup Table 1.xls"
sheet_reference <- "summary"
skip_row <- 0
coln_reference <- "patient"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- "coad"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
coad_meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))
target_recount3 <- "read"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
read_meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))
meta_recount3 <- bind_rows(coad_meta_recount3, read_meta_recount3)
coln_recount3 <- "tcga.gdc_cases.submitter_id"

# Match TCGA_ID between two tbls
# Select primary and normal samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_tumor_recount3 <- meta_recount3 %>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))
meta_normal_recount3 <- meta_recount3 %>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Solid Tissue Normal")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Remove all duplicate samples because we cannot confirm the reason
# why they selected one of the duplicate samples
lst_tcga_barcode <- meta_tumor_recount3[["tcga.gdc_cases.submitter_id"]]
dupl_barcode <- lst_tcga_barcode[duplicated(lst_tcga_barcode)]
if (length(dupl_barcode) != 0) {
  nondupl_barcode <- lst_tcga_barcode[-which(lst_tcga_barcode %in% dupl_barcode)]
} else {
  nondupl_barcode <- lst_tcga_barcode
}
meta_tumor_recount3 %<>%
  filter_at(., vars(all_of("tcga.gdc_cases.submitter_id")), all_vars(. %in% nondupl_barcode))

lst_tcga_barcode <- meta_normal_recount3[["tcga.gdc_cases.submitter_id"]]
dupl_barcode <- lst_tcga_barcode[duplicated(lst_tcga_barcode)]
if (length(dupl_barcode) != 0) {
  nondupl_barcode <- lst_tcga_barcode[-which(lst_tcga_barcode %in% dupl_barcode)]
} else {
  nondupl_barcode <- lst_tcga_barcode
}
meta_normal_recount3 %<>%
  filter_at(., vars(all_of("tcga.gdc_cases.submitter_id")), all_vars(. %in% nondupl_barcode))
meta_recount3 <- bind_rows(meta_tumor_recount3, meta_normal_recount3)

# Load COAD counts_data and tpm_data
target_recount3 <- "coad"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
coad_counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
coad_tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))

# Load READ counts_data and tpm_data
target_recount3 <- "read"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
read_counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
read_tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))

# Join two tbls
read_counts_recount3 %<>%
  .[, -which(colnames(.) == "gene_symbol")]
read_tpm_recount3 %<>%
  .[, -which(colnames(.) == "gene_symbol")]
counts_recount3 <- inner_join(coad_counts_recount3, read_counts_recount3,
  by = "ensembl_id")
tpm_recount3 <- inner_join(coad_counts_recount3, read_counts_recount3,
  by = "ensembl_id")

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Change the colnames from external_id to tcga_barcode
lst_external_id <- colnames(counts_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  counts_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

lst_tpm_external_id <- colnames(tpm_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  tpm_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data")
setwd(save_predir)
make_folder(paste0(toupper(target_reference), "_icluding_normal_tissue_nature_2012"))
setwd(paste0(toupper(target_reference), "_icluding_normal_tissue_nature_2012"))
save_dir <- file.path(save_predir, paste0(toupper(target_reference), "_icluding_normal_tissue_nature_2012"))
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- "nature11252-s2"
# Make extra folder
make_folder(cp_dir_1)
save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
files_2 <- list.files(path = path_2)
cp_files_2 <- files_2 %>%
  .[grep("\\..+$", .)]

# Copy files
setwd(save_reference_dir_2)
for (i in 1:length(cp_files_2)) {
  filename <- cp_files_2[i]
  file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
    filename))
}

# If there is an extra folder, make the folder
setwd(save_reference_dir_1)
cp_dir_1 <- "nature11252-s2 2"
# Make extra folder
make_folder(cp_dir_1)
save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
files_2 <- list.files(path = path_2)
cp_files_2 <- files_2 %>%
  .[grep("\\..+$", .)]

# Copy files
setwd(save_reference_dir_2)
for (i in 1:length(cp_files_2)) {
  filename <- cp_files_2[i]
  file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
    filename))
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    "nature11252-s2", dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "coad_read"
df_reference <- "2011-11-14592C-Sup Table 1.xls"
sheet_reference <- "summary"
skip_row <- 0
coln_reference <- "patient"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- "coad"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
coad_meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))
target_recount3 <- "read"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
read_meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))
meta_recount3 <- bind_rows(coad_meta_recount3, read_meta_recount3)
coln_recount3 <- "tcga.gdc_cases.submitter_id"

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Remove all duplicate samples because we cannot confirm the reason
# why they selected one of the duplicate samples
lst_tcga_barcode <- meta_recount3[["tcga.gdc_cases.submitter_id"]]
dupl_barcode <- lst_tcga_barcode[duplicated(lst_tcga_barcode)]
if (length(dupl_barcode) != 0) {
  nondupl_barcode <- lst_tcga_barcode[-which(lst_tcga_barcode %in% dupl_barcode)]
} else {
  nondupl_barcode <- lst_tcga_barcode
}
meta_recount3 %<>%
  filter_at(., vars(all_of("tcga.gdc_cases.submitter_id")), all_vars(. %in% nondupl_barcode))

# Load COAD counts_data and tpm_data
target_recount3 <- "coad"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
coad_counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
coad_tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))

# Load READ counts_data and tpm_data
target_recount3 <- "read"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
read_counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
read_tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))

# Join two tbls
read_counts_recount3 %<>%
  .[, -which(colnames(.) == "gene_symbol")]
read_tpm_recount3 %<>%
  .[, -which(colnames(.) == "gene_symbol")]
counts_recount3 <- inner_join(coad_counts_recount3, read_counts_recount3,
  by = "ensembl_id")
tpm_recount3 <- inner_join(coad_counts_recount3, read_counts_recount3,
  by = "ensembl_id")

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Change the colnames from external_id to tcga_barcode
lst_external_id <- colnames(counts_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  counts_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

lst_tpm_external_id <- colnames(tpm_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  tpm_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data")
setwd(save_predir)
make_folder(paste0(toupper(target_reference), "_nature_2012"))
setwd(paste0(toupper(target_reference), "_nature_2012"))
save_dir <- file.path(save_predir, paste0(toupper(target_reference), "_nature_2012"))
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- "nature11252-s2"
# Make extra folder
make_folder(cp_dir_1)
save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
files_2 <- list.files(path = path_2)
cp_files_2 <- files_2 %>%
  .[grep("\\..+$", .)]

# Copy files
setwd(save_reference_dir_2)
for (i in 1:length(cp_files_2)) {
  filename <- cp_files_2[i]
  file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
    filename))
}

# If there is an extra folder, make the folder
setwd(save_reference_dir_1)
cp_dir_1 <- "nature11252-s2 2"
# Make extra folder
make_folder(cp_dir_1)
save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
files_2 <- list.files(path = path_2)
cp_files_2 <- files_2 %>%
  .[grep("\\..+$", .)]

# Copy files
setwd(save_reference_dir_2)
for (i in 1:length(cp_files_2)) {
  filename <- cp_files_2[i]
  file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
    filename))
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

# Define recount3 parameters and load meta_data
target_recount3 <- "coad"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
coad_meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))
target_recount3 <- "read"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
read_meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))
meta_recount3 <- bind_rows(coad_meta_recount3, read_meta_recount3)
coln_recount3 <- "tcga.gdc_cases.submitter_id"

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 <- meta_recount3 %>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE"))

# Remove all duplicate samples because we cannot confirm the reason
# why they selected one of the duplicate samples
lst_tcga_barcode <- meta_recount3[["tcga.gdc_cases.submitter_id"]]
dupl_barcode <- lst_tcga_barcode[duplicated(lst_tcga_barcode)]
if (length(dupl_barcode) != 0) {
  nondupl_barcode <- lst_tcga_barcode[-which(lst_tcga_barcode %in% dupl_barcode)]
} else {
  nondupl_barcode <- lst_tcga_barcode
}
meta_recount3 %<>%
  filter_at(., vars(all_of("tcga.gdc_cases.submitter_id")), all_vars(. %in% nondupl_barcode))

# Load COAD counts_data and tpm_data
target_recount3 <- "coad"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
coad_counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
coad_tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))

# Load READ counts_data and tpm_data
target_recount3 <- "read"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
read_counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
read_tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))

# Join two tbls
read_counts_recount3 %<>%
  .[, -which(colnames(.) == "gene_symbol")]
read_tpm_recount3 %<>%
  .[, -which(colnames(.) == "gene_symbol")]
counts_recount3 <- inner_join(coad_counts_recount3, read_counts_recount3,
  by = "ensembl_id")
tpm_recount3 <- inner_join(coad_counts_recount3, read_counts_recount3,
  by = "ensembl_id")

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Change the colnames from external_id to tcga_barcode
lst_external_id <- colnames(counts_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  counts_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

lst_tpm_external_id <- colnames(tpm_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  tpm_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data")
setwd(save_predir)
make_folder(paste0(toupper("coad_read"), "_nonreference"))
setwd(paste0(toupper("coad_read"), "_nonreference"))
save_dir <- file.path(save_predir, paste0(toupper("coad_read"), "_nonreference"))
save_counts <- file.path(save_dir, paste0("coad_read", "_counts.txt"))
save_tpm <- file.path(save_dir, paste0("coad_read", "_tpm.txt"))
save_meta <- file.path(save_dir, paste0("coad_read", "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

# Define recount3 parameters and load meta_data
target_recount3 <- "dlbc"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 <- meta_recount3 %>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE"))

# Remove all duplicate samples because we cannot confirm the reason
# why they selected one of the duplicate samples
lst_tcga_barcode <- meta_recount3[["tcga.gdc_cases.submitter_id"]]
dupl_barcode <- lst_tcga_barcode[duplicated(lst_tcga_barcode)]
if (length(dupl_barcode) != 0) {
  nondupl_barcode <- lst_tcga_barcode[-which(lst_tcga_barcode %in% dupl_barcode)]
} else {
  nondupl_barcode <- lst_tcga_barcode
}
meta_recount3 %<>%
  filter_at(., vars(all_of("tcga.gdc_cases.submitter_id")), all_vars(. %in% nondupl_barcode))

# Load DLBC counts_data and tpm_data
target_recount3 <- "dlbc"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Change the colnames from external_id to tcga_barcode
lst_external_id <- colnames(counts_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  counts_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

lst_tpm_external_id <- colnames(tpm_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  tpm_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data")
setwd(save_predir)
make_folder(paste0(toupper(target_recount3), "_nonreference"))
setwd(paste0(toupper(target_recount3), "_nonreference"))
save_dir <- file.path(save_predir, paste0(toupper(target_recount3), "_nonreference"))
save_counts <- file.path(save_dir, paste0(target_recount3, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_recount3, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_recount3, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, "GBM_cell_2013",
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "gbm"
df_reference <- "table_s7.xlsx"
sheet_reference <- "Clinical Data"
skip_row <- 2
coln_reference <- "Case ID"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "tcga.gdc_cases.submitter_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Remove all duplicate samples because we cannot confirm the reason
# why they selected one of the duplicate samples
lst_tcga_barcode <- meta_recount3[["tcga.gdc_cases.submitter_id"]]
dupl_barcode <- lst_tcga_barcode[duplicated(lst_tcga_barcode)]
if (length(dupl_barcode) != 0) {
  nondupl_barcode <- lst_tcga_barcode[-which(lst_tcga_barcode %in% dupl_barcode)]
} else {
  nondupl_barcode <- lst_tcga_barcode
}
meta_recount3 %<>%
  filter_at(., vars(all_of("tcga.gdc_cases.submitter_id")), all_vars(. %in% nondupl_barcode))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Change the colnames from external_id to tcga_barcode
lst_external_id <- colnames(counts_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  counts_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

lst_tpm_external_id <- colnames(tpm_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  tpm_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data")
setwd(save_predir)
make_folder(paste0(toupper(target_reference), "_cell_2013"))
setwd(paste0(toupper(target_reference), "_cell_2013"))
save_dir <- file.path(save_predir, paste0(toupper(target_reference), "_cell_2013"))
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, "GBM_cell_2013")
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, "GBM_nature_2008",
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "gbm"
df_reference <- "supplementary_tables.xls"
sheet_reference <- "Table S1B - Individual samples"
skip_row <- 0
coln_reference <- "Case ID"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "tcga.gdc_cases.submitter_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Remove all duplicate samples because we cannot confirm the reason
# why they selected one of the duplicate samples
lst_tcga_barcode <- meta_recount3[["tcga.gdc_cases.submitter_id"]]
dupl_barcode <- lst_tcga_barcode[duplicated(lst_tcga_barcode)]
if (length(dupl_barcode) != 0) {
  nondupl_barcode <- lst_tcga_barcode[-which(lst_tcga_barcode %in% dupl_barcode)]
} else {
  nondupl_barcode <- lst_tcga_barcode
}
meta_recount3 %<>%
  filter_at(., vars(all_of("tcga.gdc_cases.submitter_id")), all_vars(. %in% nondupl_barcode))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Change the colnames from external_id to tcga_barcode
lst_external_id <- colnames(counts_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  counts_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

lst_tpm_external_id <- colnames(tpm_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  tpm_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data")
setwd(save_predir)
make_folder(paste0(toupper(target_reference), "_nature_2008"))
setwd(paste0(toupper(target_reference), "_nature_2008"))
save_dir <- file.path(save_predir, paste0(toupper(target_reference), "_nature_2008"))
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, "GBM_nature_2008")
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

# Define recount3 parameters and load meta_data
target_recount3 <- "gbm"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 <- meta_recount3 %>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE"))

# Remove all duplicate samples because we cannot confirm the reason
# why they selected one of the duplicate samples
lst_tcga_barcode <- meta_recount3[["tcga.gdc_cases.submitter_id"]]
dupl_barcode <- lst_tcga_barcode[duplicated(lst_tcga_barcode)]
if (length(dupl_barcode) != 0) {
  nondupl_barcode <- lst_tcga_barcode[-which(lst_tcga_barcode %in% dupl_barcode)]
} else {
  nondupl_barcode <- lst_tcga_barcode
}
meta_recount3 %<>%
  filter_at(., vars(all_of("tcga.gdc_cases.submitter_id")), all_vars(. %in% nondupl_barcode))

# Load GBM counts_data and tpm_data
target_recount3 <- "gbm"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Change the colnames from external_id to tcga_barcode
lst_external_id <- colnames(counts_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  counts_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

lst_tpm_external_id <- colnames(tpm_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  tpm_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data")
setwd(save_predir)
make_folder(paste0(toupper(target_recount3), "_nonreference"))
setwd(paste0(toupper(target_recount3), "_nonreference"))
save_dir <- file.path(save_predir, paste0(toupper(target_recount3), "_nonreference"))
save_counts <- file.path(save_dir, paste0(target_recount3, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_recount3, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_recount3, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    "nature14129-s2", dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "hnsc"
df_reference <- "1.1.xlsx"
sheet_reference <- "Sample_Information"
skip_row <- 0
coln_reference <- "Data Freeze Barcodes"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "tcga.gdc_cases.submitter_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Remove all duplicate samples because we cannot confirm the reason
# why they selected one of the duplicate samples
lst_tcga_barcode <- meta_recount3[["tcga.gdc_cases.submitter_id"]]
dupl_barcode <- lst_tcga_barcode[duplicated(lst_tcga_barcode)]
if (length(dupl_barcode) != 0) {
  nondupl_barcode <- lst_tcga_barcode[-which(lst_tcga_barcode %in% dupl_barcode)]
} else {
  nondupl_barcode <- lst_tcga_barcode
}
meta_recount3 %<>%
  filter_at(., vars(all_of("tcga.gdc_cases.submitter_id")), all_vars(. %in% nondupl_barcode))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Change the colnames from external_id to tcga_barcode
lst_external_id <- colnames(counts_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  counts_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

lst_tpm_external_id <- colnames(tpm_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  tpm_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "kich"
df_reference <- "table_s1.xlsx"
sheet_reference <- "by Patient"
skip_row <- 1
coln_reference <- "TCGA patient code"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "tcga.gdc_cases.submitter_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Remove all duplicate samples because we cannot confirm the reason
# why they selected one of the duplicate samples
lst_tcga_barcode <- meta_recount3[["tcga.gdc_cases.submitter_id"]]
dupl_barcode <- lst_tcga_barcode[duplicated(lst_tcga_barcode)]
if (length(dupl_barcode) != 0) {
  nondupl_barcode <- lst_tcga_barcode[-which(lst_tcga_barcode %in% dupl_barcode)]
} else {
  nondupl_barcode <- lst_tcga_barcode
}
meta_recount3 %<>%
  filter_at(., vars(all_of("tcga.gdc_cases.submitter_id")), all_vars(. %in% nondupl_barcode))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Change the colnames from external_id to tcga_barcode
lst_external_id <- colnames(counts_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  counts_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

lst_tpm_external_id <- colnames(tpm_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  tpm_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    "nature12222-s2", dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "kirc"
df_reference <- "Data_file_S1_ccRCC Freeze 1.4.1.xlsx"
sheet_reference <- "ccRCC Freeze 1.4.1"
skip_row <- 0
coln_reference <- "Participant Barcode"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "tcga.gdc_cases.submitter_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Remove all duplicate samples because we cannot confirm the reason
# why they selected one of the duplicate samples
lst_tcga_barcode <- meta_recount3[["tcga.gdc_cases.submitter_id"]]
dupl_barcode <- lst_tcga_barcode[duplicated(lst_tcga_barcode)]
if (length(dupl_barcode) != 0) {
  nondupl_barcode <- lst_tcga_barcode[-which(lst_tcga_barcode %in% dupl_barcode)]
} else {
  nondupl_barcode <- lst_tcga_barcode
}
meta_recount3 %<>%
  filter_at(., vars(all_of("tcga.gdc_cases.submitter_id")), all_vars(. %in% nondupl_barcode))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Change the colnames from external_id to tcga_barcode
lst_external_id <- colnames(counts_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  counts_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

lst_tpm_external_id <- colnames(tpm_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  tpm_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "kirp"
df_reference <- "nejmoa1505917_appendix_2.xlsx"
sheet_reference <- "KIRP Clinical Data"
skip_row <- 0
coln_reference <- "BCR Patient Barcode"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "tcga.gdc_cases.submitter_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Remove all duplicate samples because we cannot confirm the reason
# why they selected one of the duplicate samples
lst_tcga_barcode <- meta_recount3[["tcga.gdc_cases.submitter_id"]]
dupl_barcode <- lst_tcga_barcode[duplicated(lst_tcga_barcode)]
if (length(dupl_barcode) != 0) {
  nondupl_barcode <- lst_tcga_barcode[-which(lst_tcga_barcode %in% dupl_barcode)]
} else {
  nondupl_barcode <- lst_tcga_barcode
}
meta_recount3 %<>%
  filter_at(., vars(all_of("tcga.gdc_cases.submitter_id")), all_vars(. %in% nondupl_barcode))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Change the colnames from external_id to tcga_barcode
lst_external_id <- colnames(counts_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  counts_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

lst_tpm_external_id <- colnames(tpm_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  tpm_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    "laml_tcga_pub", dfname_reference)
  meta_reference <- read_tsv(reference_target_dir)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "laml"
df_reference <- "data_clinical_sample.txt"
sheet_reference <- "data_clinical_sample"
skip_row <- 4
coln_reference <- "SAMPLE_ID"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
colnames(meta_reference) <- meta_reference[4, ]
meta_reference %<>%
  .[-(1:4),]
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "sample_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))
meta_recount3 %<>%
  mutate(., sample_id = str_sub(.[["tcga.tcga_barcode"]], start = 1,
    end = 15), .before = "rail_id")

# Match TCGA_ID between two tbls
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))
meta_recount3 %<>%
  .[, -which(colnames(.) == "sample_id")]

# Remove all duplicate samples because we cannot confirm the reason
# why they selected one of the duplicate samples
lst_tcga_barcode <- na.omit(meta_recount3[["tcga.gdc_cases.submitter_id"]])
dupl_barcode <- lst_tcga_barcode[duplicated(lst_tcga_barcode)]
if (length(dupl_barcode) != 0) {
  nondupl_barcode <- lst_tcga_barcode[-which(lst_tcga_barcode %in% dupl_barcode)]
} else {
  nondupl_barcode <- meta_recount3[["tcga.gdc_cases.submitter_id"]]
}
meta_recount3 %<>%
  filter_at(., vars(all_of("tcga.gdc_cases.submitter_id")), all_vars(. %in% nondupl_barcode))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Change the colnames from external_id to tcga_barcode
lst_external_id <- colnames(counts_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  counts_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

lst_tpm_external_id <- colnames(tpm_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  tpm_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # If files_2 is including "LICENSE", "LICENSE" is include in cp_files_2
  if ("LICENSE" %in% files_2) {
    cp_files_2 %<>%
      c(., "LICENSE")
  }
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  # If files_2 is including "LICENSE", "LICENSE" is not included in cp_files_2
  if ("LICENSE" %in% files_2) {
    cp_dir_2 %<>%
      .[-which(. == "LICENSE")]
  }
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "lgg"
df_reference <- "nejmoa1402121_appendix_2.xlsx"
sheet_reference <- "Clinical.Output"
skip_row <- 0
coln_reference <- "Tumor"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "tcga.gdc_cases.submitter_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Remove all duplicate samples because we cannot confirm the reason
# why they selected one of the duplicate samples
lst_tcga_barcode <- meta_recount3[["tcga.gdc_cases.submitter_id"]]
dupl_barcode <- lst_tcga_barcode[duplicated(lst_tcga_barcode)]
if (length(dupl_barcode) != 0) {
  nondupl_barcode <- lst_tcga_barcode[-which(lst_tcga_barcode %in% dupl_barcode)]
} else {
  nondupl_barcode <- lst_tcga_barcode
}
meta_recount3 %<>%
  filter_at(., vars(all_of("tcga.gdc_cases.submitter_id")), all_vars(. %in% nondupl_barcode))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Change the colnames from external_id to tcga_barcode
lst_external_id <- colnames(counts_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  counts_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

lst_tpm_external_id <- colnames(tpm_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  tpm_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "lihc"
df_reference <- "table_s1.xlsx"
sheet_reference <- "Table S1A - core sample set"
skip_row <- 3
coln_reference <- "UUID"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "tcga.gdc_cases.samples.sample_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Remove all duplicate samples because we cannot confirm the reason
# why they selected one of the duplicate samples
lst_tcga_barcode <- meta_recount3[["tcga.gdc_cases.submitter_id"]]
dupl_barcode <- lst_tcga_barcode[duplicated(lst_tcga_barcode)]
if (length(dupl_barcode) != 0) {
  nondupl_barcode <- lst_tcga_barcode[-which(lst_tcga_barcode %in% dupl_barcode)]
} else {
  nondupl_barcode <- lst_tcga_barcode
}
meta_recount3 %<>%
  filter_at(., vars(all_of("tcga.gdc_cases.submitter_id")), all_vars(. %in% nondupl_barcode))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Change the colnames from external_id to tcga_barcode
lst_external_id <- colnames(counts_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  counts_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

lst_tpm_external_id <- colnames(tpm_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  tpm_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "luad"
df_reference <- "supplementary_tables.xlsx"
sheet_reference <- "S_Table 7-Clinical&Molec_Summar"
skip_row <- 4
coln_reference <- "Tumor ID"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "tcga.gdc_cases.submitter_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Remove all duplicate samples because we cannot confirm the reason
# why they selected one of the duplicate samples
lst_tcga_barcode <- meta_recount3[["tcga.gdc_cases.submitter_id"]]
dupl_barcode <- lst_tcga_barcode[duplicated(lst_tcga_barcode)]
if (length(dupl_barcode) != 0) {
  nondupl_barcode <- lst_tcga_barcode[-which(lst_tcga_barcode %in% dupl_barcode)]
} else {
  nondupl_barcode <- lst_tcga_barcode
}
meta_recount3 %<>%
  filter_at(., vars(all_of("tcga.gdc_cases.submitter_id")), all_vars(. %in% nondupl_barcode))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Change the colnames from external_id to tcga_barcode
lst_external_id <- colnames(counts_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  counts_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

lst_tpm_external_id <- colnames(tpm_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  tpm_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    "nature11404-s2", dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "lusc"
df_reference <- "data.file.S7.5.clinical.and.genomic.data.table.xls"
sheet_reference <- "LUSC_CpG_Filtered.patients.coun"
skip_row <- 3
coln_reference <- "Tumor ID"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- str_replace(meta_reference[[coln_reference]], pattern = "LUSC",
  replacement = "TCGA")

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "tcga.gdc_cases.submitter_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Remove all duplicate samples because we cannot confirm the reason
# why they selected one of the duplicate samples
lst_tcga_barcode <- meta_recount3[["tcga.gdc_cases.submitter_id"]]
dupl_barcode <- lst_tcga_barcode[duplicated(lst_tcga_barcode)]
if (length(dupl_barcode) != 0) {
  nondupl_barcode <- lst_tcga_barcode[-which(lst_tcga_barcode %in% dupl_barcode)]
} else {
  nondupl_barcode <- lst_tcga_barcode
}
meta_recount3 %<>%
  filter_at(., vars(all_of("tcga.gdc_cases.submitter_id")), all_vars(. %in% nondupl_barcode))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Change the colnames from external_id to tcga_barcode
lst_external_id <- colnames(counts_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  counts_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

lst_tpm_external_id <- colnames(tpm_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  tpm_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "meso"
df_reference <- "table_s1.xlsx"
sheet_reference <- "1B_MPM_Master_Patient_Table"
skip_row <- 0
coln_reference <- "TCGA_barcode"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "TCGA_barcode"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))
meta_recount3 %<>%
  mutate(., TCGA_barcode = str_sub(.[["tcga.gdc_cases.samples.submitter_id"]],
    start = 1, end = -2), .before = "rail_id")

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))
meta_recount3 %<>%
  .[, - which(colnames(.) == "TCGA_barcode")]
  
# Remove all duplicate samples because we cannot confirm the reason
# why they selected one of the duplicate samples
lst_tcga_barcode <- meta_recount3[["tcga.gdc_cases.submitter_id"]]
dupl_barcode <- lst_tcga_barcode[duplicated(lst_tcga_barcode)]
if (length(dupl_barcode) != 0) {
  nondupl_barcode <- lst_tcga_barcode[-which(lst_tcga_barcode %in% dupl_barcode)]
} else {
  nondupl_barcode <- lst_tcga_barcode
}
meta_recount3 %<>%
  filter_at(., vars(all_of("tcga.gdc_cases.submitter_id")), all_vars(. %in% nondupl_barcode))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Change the colnames from external_id to tcga_barcode
lst_external_id <- colnames(counts_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  counts_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

lst_tpm_external_id <- colnames(tpm_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  tpm_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    "nature10166-s2", dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "ov"
df_reference <- "2010-09-11380C-Table_S1.2.xlsx"
sheet_reference <- "KeyclinicalDAta"
skip_row <- 0
coln_reference <- "BCRPATIENTBARCODE"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "tcga.gdc_cases.submitter_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Remove all duplicate samples because we cannot confirm the reason
# why they selected one of the duplicate samples
lst_tcga_barcode <- meta_recount3[["tcga.gdc_cases.submitter_id"]]
dupl_barcode <- lst_tcga_barcode[duplicated(lst_tcga_barcode)]
if (length(dupl_barcode) != 0) {
  nondupl_barcode <- lst_tcga_barcode[-which(lst_tcga_barcode %in% dupl_barcode)]
} else {
  nondupl_barcode <- lst_tcga_barcode
}
meta_recount3 %<>%
  filter_at(., vars(all_of("tcga.gdc_cases.submitter_id")), all_vars(. %in% nondupl_barcode))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Change the colnames from external_id to tcga_barcode
lst_external_id <- colnames(counts_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  counts_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

lst_tpm_external_id <- colnames(tpm_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  tpm_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "paad"
df_reference <- "table_s1.xlsx"
sheet_reference <- "FreezeSamples"
skip_row <- 1
coln_reference <- "Tumor Sample ID"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "tcga.gdc_cases.samples.submitter_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Remove all duplicate samples because we cannot confirm the reason
# why they selected one of the duplicate samples
lst_tcga_barcode <- meta_recount3[["tcga.gdc_cases.submitter_id"]]
dupl_barcode <- lst_tcga_barcode[duplicated(lst_tcga_barcode)]
if (length(dupl_barcode) != 0) {
  nondupl_barcode <- lst_tcga_barcode[-which(lst_tcga_barcode %in% dupl_barcode)]
} else {
  nondupl_barcode <- lst_tcga_barcode
}
meta_recount3 %<>%
  filter_at(., vars(all_of("tcga.gdc_cases.submitter_id")), all_vars(. %in% nondupl_barcode))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Change the colnames from external_id to tcga_barcode
lst_external_id <- colnames(counts_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  counts_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

lst_tpm_external_id <- colnames(tpm_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  tpm_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "pcpg"
df_reference <- "table_s2.xls"
sheet_reference <- "Master Data"
skip_row <- 2
coln_reference <- "Sample ID"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "tcga.gdc_cases.samples.submitter_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Remove all duplicate samples because we cannot confirm the reason
# why they selected one of the duplicate samples
lst_tcga_barcode <- meta_recount3[["tcga.gdc_cases.submitter_id"]]
dupl_barcode <- lst_tcga_barcode[duplicated(lst_tcga_barcode)]
if (length(dupl_barcode) != 0) {
  nondupl_barcode <- lst_tcga_barcode[-which(lst_tcga_barcode %in% dupl_barcode)]
} else {
  nondupl_barcode <- lst_tcga_barcode
}
meta_recount3 %<>%
  filter_at(., vars(all_of("tcga.gdc_cases.submitter_id")), all_vars(. %in% nondupl_barcode))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Change the colnames from external_id to tcga_barcode
lst_external_id <- colnames(counts_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  counts_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

lst_tpm_external_id <- colnames(tpm_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  tpm_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "prad"
df_reference <- "table_s1.xls"
sheet_reference <- "Table S1A. Annotation"
skip_row <- 0
coln_reference <- "SAMPLE_ID"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "sample_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))
meta_recount3 %<>%
  mutate(., sample_id = str_sub(.[["tcga.gdc_cases.samples.submitter_id"]],
    start = 1, end = -2), .before = "rail_id")

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))
meta_recount3 %<>%
  .[, - which(colnames(.) == "sample_id")]

# Remove all duplicate samples because we cannot confirm the reason
# why they selected one of the duplicate samples
lst_tcga_barcode <- meta_recount3[["tcga.gdc_cases.submitter_id"]]
dupl_barcode <- lst_tcga_barcode[duplicated(lst_tcga_barcode)]
if (length(dupl_barcode) != 0) {
  nondupl_barcode <- lst_tcga_barcode[-which(lst_tcga_barcode %in% dupl_barcode)]
} else {
  nondupl_barcode <- lst_tcga_barcode
}
meta_recount3 %<>%
  filter_at(., vars(all_of("tcga.gdc_cases.submitter_id")), all_vars(. %in% nondupl_barcode))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Change the colnames from external_id to tcga_barcode
lst_external_id <- colnames(counts_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  counts_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

lst_tpm_external_id <- colnames(tpm_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  tpm_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "sarc_adult"
df_reference <- "table_s1.xlsx"
sheet_reference <- "Detailed sample features"
skip_row <- 1
coln_reference <- "TCGA barcode"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- "sarc"
coln_recount3 <- "TCGA_barcode"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))
meta_recount3 %<>%
  mutate(., TCGA_barcode = str_sub(.[["tcga.gdc_cases.samples.submitter_id"]],
    start = 1, end = -2), .before = "rail_id")

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))
meta_recount3 %<>%
  .[, - which(colnames(.) == "TCGA_barcode")]

# Remove all duplicate samples because we cannot confirm the reason
# why they selected one of the duplicate samples
lst_tcga_barcode <- meta_recount3[["tcga.gdc_cases.submitter_id"]]
dupl_barcode <- lst_tcga_barcode[duplicated(lst_tcga_barcode)]
if (length(dupl_barcode) != 0) {
  nondupl_barcode <- lst_tcga_barcode[-which(lst_tcga_barcode %in% dupl_barcode)]
} else {
  nondupl_barcode <- lst_tcga_barcode
}
meta_recount3 %<>%
  filter_at(., vars(all_of("tcga.gdc_cases.submitter_id")), all_vars(. %in% nondupl_barcode))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Change the colnames from external_id to tcga_barcode
lst_external_id <- colnames(counts_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  counts_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

lst_tpm_external_id <- colnames(tpm_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  tpm_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "skcm"
df_reference <- "table_s1.xlsx"
sheet_reference <- "Supplemental Table S1D"
skip_row <- 1
coln_reference <- "Name"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "Name"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))
meta_recount3 %<>%
  mutate(., Name = str_sub(.[["tcga.gdc_cases.samples.submitter_id"]],
    start = 1, end = -2), .before = "rail_id")

# Match TCGA_ID between two tbls
# Select primary  and metastatic samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor" | . == "Metastatic")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))
meta_recount3 %<>%
  .[, - which(colnames(.) == "Name")]
  
# Remove all duplicate samples because we cannot confirm the reason
# why they selected one of the duplicate samples
lst_tcga_barcode <- meta_recount3[["tcga.gdc_cases.submitter_id"]]
dupl_barcode <- lst_tcga_barcode[duplicated(lst_tcga_barcode)]
if (length(dupl_barcode) != 0) {
  nondupl_barcode <- lst_tcga_barcode[-which(lst_tcga_barcode %in% dupl_barcode)]
} else {
  nondupl_barcode <- lst_tcga_barcode
}
meta_recount3 %<>%
  filter_at(., vars(all_of("tcga.gdc_cases.submitter_id")), all_vars(. %in% nondupl_barcode))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Change the colnames from external_id to tcga_barcode
lst_external_id <- colnames(counts_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  counts_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

lst_tpm_external_id <- colnames(tpm_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  tpm_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    "stad_tcga_pub", dfname_reference)
  meta_reference <- read_tsv(reference_target_dir)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "stad"
df_reference <- "data_clinical_sample.txt"
sheet_reference <- "data_clinical_sample"
skip_row <- 4
coln_reference <- "SAMPLE_ID"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
colnames(meta_reference) <- meta_reference[4, ]
meta_reference %<>%
  .[-(1:4),]
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "sample_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))
meta_recount3 %<>%
  mutate(., sample_id = str_sub(.[["tcga.tcga_barcode"]], start = 1,
    end = 15), .before = "rail_id")

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))
meta_recount3 %<>%
  .[, -which(colnames(.) == "sample_id")]

# Remove all duplicate samples because we cannot confirm the reason
# why they selected one of the duplicate samples
lst_tcga_barcode <- na.omit(meta_recount3[["tcga.gdc_cases.submitter_id"]])
dupl_barcode <- lst_tcga_barcode[duplicated(lst_tcga_barcode)]
if (length(dupl_barcode) != 0) {
  nondupl_barcode <- lst_tcga_barcode[-which(lst_tcga_barcode %in% dupl_barcode)]
} else {
  nondupl_barcode <- meta_recount3[["tcga.gdc_cases.submitter_id"]]
}
meta_recount3 %<>%
  filter_at(., vars(all_of("tcga.gdc_cases.submitter_id")), all_vars(. %in% nondupl_barcode))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Change the colnames from external_id to tcga_barcode
lst_external_id <- colnames(counts_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  counts_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

lst_tpm_external_id <- colnames(tpm_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  tpm_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # If files_2 is including "LICENSE", "LICENSE" is include in cp_files_2
  if ("LICENSE" %in% files_2) {
    cp_files_2 %<>%
      c(., "LICENSE")
  }
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  # If files_2 is including "LICENSE", "LICENSE" is not included in cp_files_2
  if ("LICENSE" %in% files_2) {
    cp_dir_2 %<>%
      .[-which(. == "LICENSE")]
  }
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "stes"
df_reference <- "table_s1.xlsx"
sheet_reference <- "Supplementary Table 1"
skip_row <- 1
coln_reference <- "barcode"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- "esca"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
esca_meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))
target_recount3 <- "stad"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
stad_meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))
meta_recount3 <- bind_rows(esca_meta_recount3, stad_meta_recount3)
coln_recount3 <- "tcga.gdc_cases.submitter_id"

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 <- meta_recount3 %>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Remove all duplicate samples because we cannot confirm the reason
# why they selected one of the duplicate samples
lst_tcga_barcode <- meta_recount3[["tcga.gdc_cases.submitter_id"]]
dupl_barcode <- lst_tcga_barcode[duplicated(lst_tcga_barcode)]
if (length(dupl_barcode) != 0) {
  nondupl_barcode <- lst_tcga_barcode[-which(lst_tcga_barcode %in% dupl_barcode)]
} else {
  nondupl_barcode <- lst_tcga_barcode
}
meta_recount3 %<>%
  filter_at(., vars(all_of("tcga.gdc_cases.submitter_id")), all_vars(. %in% nondupl_barcode))

# Load ESCA counts_data and tpm_data
target_recount3 <- "esca"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
esca_counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
esca_tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))

# Load STAD counts_data and tpm_data
target_recount3 <- "stad"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
stad_counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
stad_tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))

# Join two tbls
stad_counts_recount3 %<>%
  .[, -which(colnames(.) == "gene_symbol")]
stad_tpm_recount3 %<>%
  .[, -which(colnames(.) == "gene_symbol")]
counts_recount3 <- inner_join(esca_counts_recount3, stad_counts_recount3,
  by = "ensembl_id")
tpm_recount3 <- inner_join(esca_counts_recount3, stad_counts_recount3,
  by = "ensembl_id")

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Change the colnames from external_id to tcga_barcode
lst_external_id <- colnames(counts_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  counts_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

lst_tpm_external_id <- colnames(tpm_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  tpm_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "tgct"
df_reference <- "table_s1.xlsx"
sheet_reference <- "TGCT_master_sif_072016"
skip_row <- 2
coln_reference <- "sample"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "sample"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))
meta_recount3 %<>%
  mutate(., sample = str_sub(.[["tcga.gdc_cases.samples.submitter_id"]],
    start = 1, end = -2), .before = "rail_id")

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))
meta_recount3 %<>%
  .[, - which(colnames(.) == "sample")]
  
# Remove all duplicate samples because we cannot confirm the reason
# why they selected one of the duplicate samples
lst_tcga_barcode <- meta_recount3[["tcga.gdc_cases.submitter_id"]]
dupl_barcode <- lst_tcga_barcode[duplicated(lst_tcga_barcode)]
if (length(dupl_barcode) != 0) {
  nondupl_barcode <- lst_tcga_barcode[-which(lst_tcga_barcode %in% dupl_barcode)]
} else {
  nondupl_barcode <- lst_tcga_barcode
}
meta_recount3 %<>%
  filter_at(., vars(all_of("tcga.gdc_cases.submitter_id")), all_vars(. %in% nondupl_barcode))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Change the colnames from external_id to tcga_barcode
lst_external_id <- colnames(counts_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  counts_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

lst_tpm_external_id <- colnames(tpm_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  tpm_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "thca_ptc"
df_reference <- "table_s2.xlsx"
sheet_reference <- "THCA-TP (496) "
skip_row <- 0
coln_reference <- "sample"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- "thca"
coln_recount3 <- "tcga.gdc_cases.samples.submitter_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Remove all duplicate samples because we cannot confirm the reason
# why they selected one of the duplicate samples
lst_tcga_barcode <- meta_recount3[["tcga.gdc_cases.submitter_id"]]
dupl_barcode <- lst_tcga_barcode[duplicated(lst_tcga_barcode)]
if (length(dupl_barcode) != 0) {
  nondupl_barcode <- lst_tcga_barcode[-which(lst_tcga_barcode %in% dupl_barcode)]
} else {
  nondupl_barcode <- lst_tcga_barcode
}
meta_recount3 %<>%
  filter_at(., vars(all_of("tcga.gdc_cases.submitter_id")), all_vars(. %in% nondupl_barcode))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Change the colnames from external_id to tcga_barcode
lst_external_id <- colnames(counts_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  counts_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

lst_tpm_external_id <- colnames(tpm_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  tpm_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "thym_tets"
df_reference <- "table_s1.xlsx"
sheet_reference <- "Filtered Freeze Set"
skip_row <- 1
coln_reference <- "UUID"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- "thym"
coln_recount3 <- "tcga.gdc_cases.samples.sample_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Remove all duplicate samples because we cannot confirm the reason
# why they selected one of the duplicate samples
lst_tcga_barcode <- meta_recount3[["tcga.gdc_cases.submitter_id"]]
dupl_barcode <- lst_tcga_barcode[duplicated(lst_tcga_barcode)]
if (length(dupl_barcode) != 0) {
  nondupl_barcode <- lst_tcga_barcode[-which(lst_tcga_barcode %in% dupl_barcode)]
} else {
  nondupl_barcode <- lst_tcga_barcode
}
meta_recount3 %<>%
  filter_at(., vars(all_of("tcga.gdc_cases.submitter_id")), all_vars(. %in% nondupl_barcode))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Change the colnames from external_id to tcga_barcode
lst_external_id <- colnames(counts_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  counts_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

lst_tpm_external_id <- colnames(tpm_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  tpm_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    "supplementary_data", dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "ucec"
df_reference <- "datafile.S1.1.KeyClinicalData.xls"
sheet_reference <- "373cases"
skip_row <- 0
coln_reference <- "bcr_patient_barcode"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "tcga.gdc_cases.submitter_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Remove all duplicate samples because we cannot confirm the reason
# why they selected one of the duplicate samples
lst_tcga_barcode <- meta_recount3[["tcga.gdc_cases.submitter_id"]]
dupl_barcode <- lst_tcga_barcode[duplicated(lst_tcga_barcode)]
if (length(dupl_barcode) != 0) {
  nondupl_barcode <- lst_tcga_barcode[-which(lst_tcga_barcode %in% dupl_barcode)]
} else {
  nondupl_barcode <- lst_tcga_barcode
}
meta_recount3 %<>%
  filter_at(., vars(all_of("tcga.gdc_cases.submitter_id")), all_vars(. %in% nondupl_barcode))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Change the colnames from external_id to tcga_barcode
lst_external_id <- colnames(counts_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  counts_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

lst_tpm_external_id <- colnames(tpm_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  tpm_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "ucs"
df_reference <- "table_s1.xlsx"
sheet_reference <- "TableS1"
skip_row <- 1
coln_reference <- "bcr_patient_barcode"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "tcga.gdc_cases.submitter_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Remove all duplicate samples because we cannot confirm the reason
# why they selected one of the duplicate samples
lst_tcga_barcode <- meta_recount3[["tcga.gdc_cases.submitter_id"]]
dupl_barcode <- lst_tcga_barcode[duplicated(lst_tcga_barcode)]
if (length(dupl_barcode) != 0) {
  nondupl_barcode <- lst_tcga_barcode[-which(lst_tcga_barcode %in% dupl_barcode)]
} else {
  nondupl_barcode <- lst_tcga_barcode
}
meta_recount3 %<>%
  filter_at(., vars(all_of("tcga.gdc_cases.submitter_id")), all_vars(. %in% nondupl_barcode))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Change the colnames from external_id to tcga_barcode
lst_external_id <- colnames(counts_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  counts_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

lst_tpm_external_id <- colnames(tpm_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  tpm_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "uvm"
df_reference <- "table_s1.xlsx"
sheet_reference <- "Summary_by_Case"
skip_row <- 3
coln_reference <- "Patient ID"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "tcga.gdc_cases.submitter_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Remove all duplicate samples because we cannot confirm the reason
# why they selected one of the duplicate samples
lst_tcga_barcode <- meta_recount3[["tcga.gdc_cases.submitter_id"]]
dupl_barcode <- lst_tcga_barcode[duplicated(lst_tcga_barcode)]
if (length(dupl_barcode) != 0) {
  nondupl_barcode <- lst_tcga_barcode[-which(lst_tcga_barcode %in% dupl_barcode)]
} else {
  nondupl_barcode <- lst_tcga_barcode
}
meta_recount3 %<>%
  filter_at(., vars(all_of("tcga.gdc_cases.submitter_id")), all_vars(. %in% nondupl_barcode))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Change the colnames from external_id to tcga_barcode
lst_external_id <- colnames(counts_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  counts_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

lst_tpm_external_id <- colnames(tpm_recount3) %>%
  .[-c(which(. == "ensembl_id"), which(. == "gene_symbol"))]
for (i in 1:length(lst_external_id)) {
  tcga_barcode <- filter_at(meta_recount3, vars(all_of("external_id")),
    all_vars(. == lst_external_id[i]))[["tcga.tcga_barcode"]]
  tpm_recount3 %<>%
    rename_at(., vars(all_of(lst_external_id[i])), ~tcga_barcode)
}

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Load package
library(tidyverse)
library(magrittr)
library(readxl)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "acc"
df_reference <- "table_s1.xlsx"
sheet_reference <- "Data Overview"
skip_row <- 6
coln_reference <- "TCGA_ID"
except_row <- "Platform"
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "tcga.gdc_cases.submitter_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data_including_duplicated")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "blca_mibc"
df_reference <- "table_s1.xlsx"
sheet_reference <- "Master table"
skip_row <- 0
coln_reference <- "BCR patient uuid"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- tolower(meta_reference[[coln_reference]])

# Define recount3 parameters and load meta_data
target_recount3 <- "blca"
coln_recount3 <- "tcga.gdc_cases.case_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data_including_duplicated")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "brca"
df_reference <- "table_s1.xlsx"
sheet_reference <- "Suppl. Table 1"
skip_row <- 2
coln_reference <- "mRNA"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "tcga.tcga_barcode"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data_including_duplicated")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "chol"
df_reference <- "tables.xlsx"
sheet_reference <- "ST1. CHOL Sample Info"
skip_row <- 3
coln_reference <- "...1"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
meta_reference %<>% filter_at(., vars(all_of("Sample Type")), all_vars(. == "Primary solid Tumor"))
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "tcga.gdc_cases.submitter_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data_including_duplicated")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    "nature11252-s2", dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "coad_read"
df_reference <- "2011-11-14592C-Sup Table 1.xls"
sheet_reference <- "summary"
skip_row <- 0
coln_reference <- "patient"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- "coad"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
coad_meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))
target_recount3 <- "read"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
read_meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))
meta_recount3 <- bind_rows(coad_meta_recount3, read_meta_recount3)
coln_recount3 <- "tcga.gdc_cases.submitter_id"

# Match TCGA_ID between two tbls
# Select primary and normal samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_tumor_recount3 <- meta_recount3 %>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))
meta_normal_recount3 <- meta_recount3 %>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Solid Tissue Normal")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Load COAD counts_data and tpm_data
target_recount3 <- "coad"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
coad_counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
coad_tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))

# Load READ counts_data and tpm_data
target_recount3 <- "read"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
read_counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
read_tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))

# Join two tbls
read_counts_recount3 %<>%
  .[, -which(colnames(.) == "gene_symbol")]
read_tpm_recount3 %<>%
  .[, -which(colnames(.) == "gene_symbol")]
counts_recount3 <- inner_join(coad_counts_recount3, read_counts_recount3,
  by = "ensembl_id")
tpm_recount3 <- inner_join(coad_counts_recount3, read_counts_recount3,
  by = "ensembl_id")

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data_including_duplicated")
setwd(save_predir)
make_folder(paste0(toupper(target_reference), "_icluding_normal_tissue_nature_2012"))
setwd(paste0(toupper(target_reference), "_icluding_normal_tissue_nature_2012"))
save_dir <- file.path(save_predir, paste0(toupper(target_reference), "_icluding_normal_tissue_nature_2012"))
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- "nature11252-s2"
# Make extra folder
make_folder(cp_dir_1)
save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
files_2 <- list.files(path = path_2)
cp_files_2 <- files_2 %>%
  .[grep("\\..+$", .)]

# Copy files
setwd(save_reference_dir_2)
for (i in 1:length(cp_files_2)) {
  filename <- cp_files_2[i]
  file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
    filename))
}

# If there is an extra folder, make the folder
setwd(save_reference_dir_1)
cp_dir_1 <- "nature11252-s2 2"
# Make extra folder
make_folder(cp_dir_1)
save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
files_2 <- list.files(path = path_2)
cp_files_2 <- files_2 %>%
  .[grep("\\..+$", .)]

# Copy files
setwd(save_reference_dir_2)
for (i in 1:length(cp_files_2)) {
  filename <- cp_files_2[i]
  file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
    filename))
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    "nature11252-s2", dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "coad_read"
df_reference <- "2011-11-14592C-Sup Table 1.xls"
sheet_reference <- "summary"
skip_row <- 0
coln_reference <- "patient"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- "coad"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
coad_meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))
target_recount3 <- "read"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
read_meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))
meta_recount3 <- bind_rows(coad_meta_recount3, read_meta_recount3)
coln_recount3 <- "tcga.gdc_cases.submitter_id"

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Load COAD counts_data and tpm_data
target_recount3 <- "coad"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
coad_counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
coad_tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))

# Load READ counts_data and tpm_data
target_recount3 <- "read"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
read_counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
read_tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))

# Join two tbls
read_counts_recount3 %<>%
  .[, -which(colnames(.) == "gene_symbol")]
read_tpm_recount3 %<>%
  .[, -which(colnames(.) == "gene_symbol")]
counts_recount3 <- inner_join(coad_counts_recount3, read_counts_recount3,
  by = "ensembl_id")
tpm_recount3 <- inner_join(coad_counts_recount3, read_counts_recount3,
  by = "ensembl_id")

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data_including_duplicated")
setwd(save_predir)
make_folder(paste0(toupper(target_reference), "_nature_2012"))
setwd(paste0(toupper(target_reference), "_nature_2012"))
save_dir <- file.path(save_predir, paste0(toupper(target_reference), "_nature_2012"))
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- "nature11252-s2"
# Make extra folder
make_folder(cp_dir_1)
save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
files_2 <- list.files(path = path_2)
cp_files_2 <- files_2 %>%
  .[grep("\\..+$", .)]

# Copy files
setwd(save_reference_dir_2)
for (i in 1:length(cp_files_2)) {
  filename <- cp_files_2[i]
  file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
    filename))
}

# If there is an extra folder, make the folder
setwd(save_reference_dir_1)
cp_dir_1 <- "nature11252-s2 2"
# Make extra folder
make_folder(cp_dir_1)
save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
files_2 <- list.files(path = path_2)
cp_files_2 <- files_2 %>%
  .[grep("\\..+$", .)]

# Copy files
setwd(save_reference_dir_2)
for (i in 1:length(cp_files_2)) {
  filename <- cp_files_2[i]
  file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
    filename))
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

# Define recount3 parameters and load meta_data
target_recount3 <- "coad"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
coad_meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))
target_recount3 <- "read"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
read_meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))
meta_recount3 <- bind_rows(coad_meta_recount3, read_meta_recount3)
coln_recount3 <- "tcga.gdc_cases.submitter_id"

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 <- meta_recount3 %>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE"))

# Load COAD counts_data and tpm_data
target_recount3 <- "coad"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
coad_counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
coad_tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))

# Load READ counts_data and tpm_data
target_recount3 <- "read"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
read_counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
read_tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))

# Join two tbls
read_counts_recount3 %<>%
  .[, -which(colnames(.) == "gene_symbol")]
read_tpm_recount3 %<>%
  .[, -which(colnames(.) == "gene_symbol")]
counts_recount3 <- inner_join(coad_counts_recount3, read_counts_recount3,
  by = "ensembl_id")
tpm_recount3 <- inner_join(coad_counts_recount3, read_counts_recount3,
  by = "ensembl_id")

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data_including_duplicated")
setwd(save_predir)
make_folder(paste0(toupper("coad_read"), "_nonreference"))
setwd(paste0(toupper("coad_read"), "_nonreference"))
save_dir <- file.path(save_predir, paste0(toupper("coad_read"), "_nonreference"))
save_counts <- file.path(save_dir, paste0("coad_read", "_counts.txt"))
save_tpm <- file.path(save_dir, paste0("coad_read", "_tpm.txt"))
save_meta <- file.path(save_dir, paste0("coad_read", "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

# Define recount3 parameters and load meta_data
target_recount3 <- "dlbc"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 <- meta_recount3 %>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE"))

# Load DLBC counts_data and tpm_data
target_recount3 <- "dlbc"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data_including_duplicated")
setwd(save_predir)
make_folder(paste0(toupper(target_recount3), "_nonreference"))
setwd(paste0(toupper(target_recount3), "_nonreference"))
save_dir <- file.path(save_predir, paste0(toupper(target_recount3), "_nonreference"))
save_counts <- file.path(save_dir, paste0(target_recount3, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_recount3, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_recount3, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, "GBM_cell_2013",
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "gbm"
df_reference <- "table_s7.xlsx"
sheet_reference <- "Clinical Data"
skip_row <- 2
coln_reference <- "Case ID"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "tcga.gdc_cases.submitter_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data_including_duplicated")
setwd(save_predir)
make_folder(paste0(toupper(target_reference), "_cell_2013"))
setwd(paste0(toupper(target_reference), "_cell_2013"))
save_dir <- file.path(save_predir, paste0(toupper(target_reference), "_cell_2013"))
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, "GBM_cell_2013")
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, "GBM_nature_2008",
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "gbm"
df_reference <- "supplementary_tables.xls"
sheet_reference <- "Table S1B - Individual samples"
skip_row <- 0
coln_reference <- "Case ID"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "tcga.gdc_cases.submitter_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data_including_duplicated")
setwd(save_predir)
make_folder(paste0(toupper(target_reference), "_nature_2008"))
setwd(paste0(toupper(target_reference), "_nature_2008"))
save_dir <- file.path(save_predir, paste0(toupper(target_reference), "_nature_2008"))
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, "GBM_nature_2008")
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

# Define recount3 parameters and load meta_data
target_recount3 <- "gbm"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 <- meta_recount3 %>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE"))

# Load GBM counts_data and tpm_data
target_recount3 <- "gbm"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data_including_duplicated")
setwd(save_predir)
make_folder(paste0(toupper(target_recount3), "_nonreference"))
setwd(paste0(toupper(target_recount3), "_nonreference"))
save_dir <- file.path(save_predir, paste0(toupper(target_recount3), "_nonreference"))
save_counts <- file.path(save_dir, paste0(target_recount3, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_recount3, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_recount3, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    "nature14129-s2", dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "hnsc"
df_reference <- "1.1.xlsx"
sheet_reference <- "Sample_Information"
skip_row <- 0
coln_reference <- "Data Freeze Barcodes"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "tcga.gdc_cases.submitter_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data_including_duplicated")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "kich"
df_reference <- "table_s1.xlsx"
sheet_reference <- "by Patient"
skip_row <- 1
coln_reference <- "TCGA patient code"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "tcga.gdc_cases.submitter_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data_including_duplicated")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    "nature12222-s2", dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "kirc"
df_reference <- "Data_file_S1_ccRCC Freeze 1.4.1.xlsx"
sheet_reference <- "ccRCC Freeze 1.4.1"
skip_row <- 0
coln_reference <- "Participant Barcode"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "tcga.gdc_cases.submitter_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data_including_duplicated")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "kirp"
df_reference <- "nejmoa1505917_appendix_2.xlsx"
sheet_reference <- "KIRP Clinical Data"
skip_row <- 0
coln_reference <- "BCR Patient Barcode"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "tcga.gdc_cases.submitter_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data_including_duplicated")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    "laml_tcga_pub", dfname_reference)
  meta_reference <- read_tsv(reference_target_dir)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "laml"
df_reference <- "data_clinical_sample.txt"
sheet_reference <- "data_clinical_sample"
skip_row <- 4
coln_reference <- "SAMPLE_ID"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
colnames(meta_reference) <- meta_reference[4, ]
meta_reference %<>%
  .[-(1:4),]
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "sample_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))
meta_recount3 %<>%
  mutate(., sample_id = str_sub(.[["tcga.tcga_barcode"]], start = 1,
    end = 15), .before = "rail_id")

# Match TCGA_ID between two tbls
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))
meta_recount3 %<>%
  .[, -which(colnames(.) == "sample_id")]

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data_including_duplicated")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # If files_2 is including "LICENSE", "LICENSE" is include in cp_files_2
  if ("LICENSE" %in% files_2) {
    cp_files_2 %<>%
      c(., "LICENSE")
  }
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  # If files_2 is including "LICENSE", "LICENSE" is not included in cp_files_2
  if ("LICENSE" %in% files_2) {
    cp_dir_2 %<>%
      .[-which(. == "LICENSE")]
  }
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "lgg"
df_reference <- "nejmoa1402121_appendix_2.xlsx"
sheet_reference <- "Clinical.Output"
skip_row <- 0
coln_reference <- "Tumor"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "tcga.gdc_cases.submitter_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data_including_duplicated")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "lihc"
df_reference <- "table_s1.xlsx"
sheet_reference <- "Table S1A - core sample set"
skip_row <- 3
coln_reference <- "UUID"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "tcga.gdc_cases.samples.sample_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data_including_duplicated")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "luad"
df_reference <- "supplementary_tables.xlsx"
sheet_reference <- "S_Table 7-Clinical&Molec_Summar"
skip_row <- 4
coln_reference <- "Tumor ID"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "tcga.gdc_cases.submitter_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data_including_duplicated")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    "nature11404-s2", dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "lusc"
df_reference <- "data.file.S7.5.clinical.and.genomic.data.table.xls"
sheet_reference <- "LUSC_CpG_Filtered.patients.coun"
skip_row <- 3
coln_reference <- "Tumor ID"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- str_replace(meta_reference[[coln_reference]], pattern = "LUSC",
  replacement = "TCGA")

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "tcga.gdc_cases.submitter_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data_including_duplicated")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "meso"
df_reference <- "table_s1.xlsx"
sheet_reference <- "1B_MPM_Master_Patient_Table"
skip_row <- 0
coln_reference <- "TCGA_barcode"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "TCGA_barcode"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))
meta_recount3 %<>%
  mutate(., TCGA_barcode = str_sub(.[["tcga.gdc_cases.samples.submitter_id"]],
    start = 1, end = -2), .before = "rail_id")

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))
meta_recount3 %<>%
  .[, - which(colnames(.) == "TCGA_barcode")]

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data_including_duplicated")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    "nature10166-s2", dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "ov"
df_reference <- "2010-09-11380C-Table_S1.2.xlsx"
sheet_reference <- "KeyclinicalDAta"
skip_row <- 0
coln_reference <- "BCRPATIENTBARCODE"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "tcga.gdc_cases.submitter_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data_including_duplicated")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "paad"
df_reference <- "table_s1.xlsx"
sheet_reference <- "FreezeSamples"
skip_row <- 1
coln_reference <- "Tumor Sample ID"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "tcga.gdc_cases.samples.submitter_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data_including_duplicated")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "pcpg"
df_reference <- "table_s2.xls"
sheet_reference <- "Master Data"
skip_row <- 2
coln_reference <- "Sample ID"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "tcga.gdc_cases.samples.submitter_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data_including_duplicated")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "prad"
df_reference <- "table_s1.xls"
sheet_reference <- "Table S1A. Annotation"
skip_row <- 0
coln_reference <- "SAMPLE_ID"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "sample_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))
meta_recount3 %<>%
  mutate(., sample_id = str_sub(.[["tcga.gdc_cases.samples.submitter_id"]],
    start = 1, end = -2), .before = "rail_id")

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))
meta_recount3 %<>%
  .[, - which(colnames(.) == "sample_id")]

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data_including_duplicated")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "sarc_adult"
df_reference <- "table_s1.xlsx"
sheet_reference <- "Detailed sample features"
skip_row <- 1
coln_reference <- "TCGA barcode"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- "sarc"
coln_recount3 <- "TCGA_barcode"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))
meta_recount3 %<>%
  mutate(., TCGA_barcode = str_sub(.[["tcga.gdc_cases.samples.submitter_id"]],
    start = 1, end = -2), .before = "rail_id")

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))
meta_recount3 %<>%
  .[, - which(colnames(.) == "TCGA_barcode")]

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data_including_duplicated")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "skcm"
df_reference <- "table_s1.xlsx"
sheet_reference <- "Supplemental Table S1D"
skip_row <- 1
coln_reference <- "Name"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "Name"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))
meta_recount3 %<>%
  mutate(., Name = str_sub(.[["tcga.gdc_cases.samples.submitter_id"]],
    start = 1, end = -2), .before = "rail_id")

# Match TCGA_ID between two tbls
# Select primary  and metastatic samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor" | . == "Metastatic")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))
meta_recount3 %<>%
  .[, - which(colnames(.) == "Name")]

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data_including_duplicated")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    "stad_tcga_pub", dfname_reference)
  meta_reference <- read_tsv(reference_target_dir)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "stad"
df_reference <- "data_clinical_sample.txt"
sheet_reference <- "data_clinical_sample"
skip_row <- 4
coln_reference <- "SAMPLE_ID"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
colnames(meta_reference) <- meta_reference[4, ]
meta_reference %<>%
  .[-(1:4),]
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "sample_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))
meta_recount3 %<>%
  mutate(., sample_id = str_sub(.[["tcga.tcga_barcode"]], start = 1,
    end = 15), .before = "rail_id")

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))
meta_recount3 %<>%
  .[, -which(colnames(.) == "sample_id")]

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data_including_duplicated")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # If files_2 is including "LICENSE", "LICENSE" is include in cp_files_2
  if ("LICENSE" %in% files_2) {
    cp_files_2 %<>%
      c(., "LICENSE")
  }
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  # If files_2 is including "LICENSE", "LICENSE" is not included in cp_files_2
  if ("LICENSE" %in% files_2) {
    cp_dir_2 %<>%
      .[-which(. == "LICENSE")]
  }
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "stes"
df_reference <- "table_s1.xlsx"
sheet_reference <- "Supplementary Table 1"
skip_row <- 1
coln_reference <- "barcode"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- "esca"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
esca_meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))
target_recount3 <- "stad"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
stad_meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))
meta_recount3 <- bind_rows(esca_meta_recount3, stad_meta_recount3)
coln_recount3 <- "tcga.gdc_cases.submitter_id"

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 <- meta_recount3 %>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Load ESCA counts_data and tpm_data
target_recount3 <- "esca"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
esca_counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
esca_tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))

# Load STAD counts_data and tpm_data
target_recount3 <- "stad"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
stad_counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
stad_tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))

# Join two tbls
stad_counts_recount3 %<>%
  .[, -which(colnames(.) == "gene_symbol")]
stad_tpm_recount3 %<>%
  .[, -which(colnames(.) == "gene_symbol")]
counts_recount3 <- inner_join(esca_counts_recount3, stad_counts_recount3,
  by = "ensembl_id")
tpm_recount3 <- inner_join(esca_counts_recount3, stad_counts_recount3,
  by = "ensembl_id")

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data_including_duplicated")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "tgct"
df_reference <- "table_s1.xlsx"
sheet_reference <- "TGCT_master_sif_072016"
skip_row <- 2
coln_reference <- "sample"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "sample"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))
meta_recount3 %<>%
  mutate(., sample = str_sub(.[["tcga.gdc_cases.samples.submitter_id"]],
    start = 1, end = -2), .before = "rail_id")

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))
meta_recount3 %<>%
  .[, - which(colnames(.) == "sample")]

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data_including_duplicated")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "thca_ptc"
df_reference <- "table_s2.xlsx"
sheet_reference <- "THCA-TP (496) "
skip_row <- 0
coln_reference <- "sample"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- "thca"
coln_recount3 <- "tcga.gdc_cases.samples.submitter_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data_including_duplicated")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "thym_tets"
df_reference <- "table_s1.xlsx"
sheet_reference <- "Filtered Freeze Set"
skip_row <- 1
coln_reference <- "UUID"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- "thym"
coln_recount3 <- "tcga.gdc_cases.samples.sample_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data_including_duplicated")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    "supplementary_data", dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "ucec"
df_reference <- "datafile.S1.1.KeyClinicalData.xls"
sheet_reference <- "373cases"
skip_row <- 0
coln_reference <- "bcr_patient_barcode"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "tcga.gdc_cases.submitter_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data_including_duplicated")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "ucs"
df_reference <- "table_s1.xlsx"
sheet_reference <- "TableS1"
skip_row <- 1
coln_reference <- "bcr_patient_barcode"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "tcga.gdc_cases.submitter_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data_including_duplicated")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

# Define directory
reference_dir <- "/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Reference_Data"
recount3_dir <- "/Volumes/G_DRIVEmobile/recount3_Rdata/TCGA"

# Functions
make_folder <- function(folder_name) {
  if (file.exists(folder_name) == FALSE) {
    dir.create(folder_name)
  }
}

make_meta_reference <- function(dfname_reference) {
  reference_target_dir <- file.path(reference_dir, toupper(target_reference),
    dfname_reference)
  meta_reference <- read_excel(reference_target_dir, sheet = sheet_reference,
    skip = skip_row)
  if (except_row != " ") {
    meta_reference %<>%
      filter_at(., vars(all_of(coln_reference)), all_vars(. != except_row))
  } else {
    meta_reference <- meta_reference
  }
}

# Define reference parameters and load meta_data
target_reference <- "uvm"
df_reference <- "table_s1.xlsx"
sheet_reference <- "Summary_by_Case"
skip_row <- 3
coln_reference <- "Patient ID"
except_row <- " "
meta_reference <- make_meta_reference(df_reference)
lst_reference_id <- meta_reference[[coln_reference]]

# Define recount3 parameters and load meta_data
target_recount3 <- target_reference
coln_recount3 <- "tcga.gdc_cases.submitter_id"
recount3_target_dir <- file.path(recount3_dir, toupper(target_recount3))
meta_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_meta.txt")))

# Match TCGA_ID between two tbls
# Select primary samples
# Remove ffpe samples
primary_recount3 <- "tcga.gdc_cases.samples.sample_type"
ffpe_recount3 <- "tcga.gdc_cases.samples.is_ffpe"
meta_recount3 %<>%
  filter_at(., vars(all_of(primary_recount3)), all_vars(. == "Primary Tumor")) %>%
  filter_at(., vars(all_of(ffpe_recount3)), all_vars(. == "FALSE")) %>%
  filter_at(., vars(all_of(coln_recount3)), all_vars(. %in% lst_reference_id))

# Load counts_data and tpm_data
counts_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_counts.txt")))
tpm_recount3 <- read_tsv(file.path(recount3_target_dir, paste0(target_recount3,
  "_tpm.txt")))

# Match TCGA_ID between counts_data and tpm_data and meta_data
counts_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))
tpm_recount3 %<>%
  select(., c(ensembl_id, gene_symbol, meta_recount3[["external_id"]]))

# Save files
save_predir <- file.path("/Volumes/G_DRIVEmobile/Revised_recount3_Rdata/Curated_Data_including_duplicated")
setwd(save_predir)
make_folder(toupper(target_reference))
save_dir <- file.path(save_predir, toupper(target_reference))
setwd(save_dir)
save_counts <- file.path(save_dir, paste0(target_reference, "_counts.txt"))
save_tpm <- file.path(save_dir, paste0(target_reference, "_tpm.txt"))
save_meta <- file.path(save_dir, paste0(target_reference, "_meta.txt"))
write_tsv(counts_recount3, save_counts)
write_tsv(tpm_recount3, save_tpm)
write_tsv(meta_recount3, save_meta)

# Copy reference
make_folder("reference")
save_reference_dir_1 <- file.path(save_dir, "reference")
setwd(save_reference_dir_1)
path_1 <- file.path(reference_dir, toupper(target_reference))
files_1 <- list.files(path = path_1)
cp_files_1 <- files_1[grep("\\..+$", files_1)]

# Exclude zip and tar files
if (length(cp_files_1[grep("\\.zip$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.zip$", .)]
}
if (length(cp_files_1[grep("\\.tar$", cp_files_1)]) != 0) {
  cp_files_1 %<>%
    .[-grep("\\.tar$", .)]
}

# Copy files
for (i in 1:length(cp_files_1)) {
  filename <- cp_files_1[i]
  file.copy(from = file.path(path_1, filename), to = file.path(save_reference_dir_1,
    filename))
}

# If there is an extra folder, make the folder
cp_dir_1 <- files_1[-grep("\\..+$", files_1)]
if (length(cp_dir_1) != 0) {
  # Make extra folder
  make_folder(cp_dir_1)
  save_reference_dir_2 <- file.path(save_reference_dir_1, cp_dir_1)
  path_2 <- file.path(reference_dir, toupper(target_reference), cp_dir_1)
  files_2 <- list.files(path = path_2)
  cp_files_2 <- files_2 %>%
    .[grep("\\..+$", .)]
  # Exclude zip and tar files
  if (length(cp_files_2[grep("\\.zip$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.zip$", .)]
  }
  if (length(cp_files_2[grep("\\.tar$", cp_files_2)]) != 0) {
    cp_files_2 %<>%
      .[-grep("\\.tar$", .)]
  }
  # Copy files
  setwd(save_reference_dir_2)
  for (i in 1:length(cp_files_2)) {
    filename <- cp_files_2[i]
    file.copy(from = file.path(path_2, filename), to = file.path(save_reference_dir_2,
      filename))
  }
  # If there is an extra folder, make the folder
  cp_dir_2 <- files_2[-grep("\\..+$", files_2)]
  if (length(cp_dir_2) != 0) {
    # Make extra folder
    make_folder(cp_dir_2)
    save_reference_dir_3 <- file.path(save_reference_dir_2, cp_dir_2)
    path_3 <- file.path(reference_dir, toupper(target_reference), cp_dir_1,
      cp_dir_2)
    files_3 <- list.files(path = path_3)
    cp_files_3 <- files_3 %>%
      .[grep("\\..+$", .)]
    # Exclude zip and tar files
    if (length(cp_files_3[grep("\\.zip$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.zip$", .)]
    }
    if (length(cp_files_3[grep("\\.tar$", cp_files_3)]) != 0) {
      cp_files_3 %<>%
        .[-grep("\\.tar$", .)]
    }
    # Copy files
    setwd(save_reference_dir_3)
    for (i in 1:length(cp_files_3)) {
      filename <- cp_files_3[i]
      file.copy(from = file.path(path_3, filename), to = file.path(save_reference_dir_3,
        filename))
    }
  }
}

traceback()

