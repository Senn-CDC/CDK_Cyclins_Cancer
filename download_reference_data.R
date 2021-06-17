# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

setwd("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/reference_data")
file_exists <- file.exists("ACC")
if (file_exists == FALSE) {
  dir.create("ACC")
}
setwd("ACC")

# Download files
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S153561081630160X-mmc1.pdf"
download.file(curl, "document_s1.pdf")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S153561081630160X-mmc2.xlsx"
download.file(curl, "table_s1.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S153561081630160X-mmc3.xlsx"
download.file(curl, "table_s2.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S153561081630160X-mmc4.xlsx"
download.file(curl, "table_s3.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S153561081630160X-mmc5.xlsx"
download.file(curl, "table_s4.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S153561081630160X-mmc6.pdf"
download.file(curl, "Comprehensive Pan-Genomic Characterization of Adrenocortical Carcinoma.pdf")

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

setwd("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/reference_data")
file_exists <- file.exists("BLCA_MIBC")
if (file_exists == FALSE) {
  dir.create("BLCA_MIBC")
}
setwd("BLCA_MIBC")

# Download files
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867417310565-mmc1.xlsx"
download.file(curl, "table_s1.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867417310565-mmc2.xlsx"
download.file(curl, "table_s2.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867417310565-mmc3.xlsx"
download.file(curl, "table_s3.xlsx")
curl <- "https://www.sciencedirect.com/science/article/pii/S0092867417310565/pdfft?isDTMRedir=true&download=true"
download.file(curl, "Comprehensive Molecular Characterization of Muscle-Invasive Bladder Cancer.pdf")

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

setwd("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/reference_data")
file_exists <- file.exists("BRCA")
if (file_exists == FALSE) {
  dir.create("BRCA")
}
setwd("BRCA")

# Download files
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867415011952-mmc1.pdf"
download.file(curl, "document_s1.pdf")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867415011952-mmc2.xlsx"
download.file(curl, "table_s1.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867415011952-mmc3.xlsx"
download.file(curl, "table_s2.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867415011952-mmc4.xlsx"
download.file(curl, "table_s3.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867415011952-mmc5.xlsx"
download.file(curl, "table_s4.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867415011952-mmc6.xlsx"
download.file(curl, "table_s5.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867415011952-mmc7.xlsx"
download.file(curl, "table_s6.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867415011952-mmc8.xlsx"
download.file(curl, "table_s7.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867415011952-mmc9.xlsx"
download.file(curl, "table_s8.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867415011952-mmc10.xlsx"
download.file(curl, "table_s9.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867415011952-mmc11.pdf"
download.file(curl, "Comprehensive Molecular Portraits of Invasive Lobular Breast Cancer.pdf")

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

setwd("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/reference_data")
file_exists <- file.exists("CESC")
if (file_exists == FALSE) {
  dir.create("CESC")
}
setwd("CESC")

# Download files
curl <- "https://static-content.springer.com/esm/art%3A10.1038%2Fnature21386/MediaObjects/41586_2017_BFnature21386_MOESM213_ESM.zip"
download.file(curl, "supplementary_data.zip")
curl <- "https://www.nature.com/articles/nature21386.pdf"
download.file(curl, "Integrated genomic and molecular characterization of cervical cancer.pdf")

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

setwd("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/reference_data")
file_exists <- file.exists("CHOL")
if (file_exists == FALSE) {
  dir.create("CHOL")
}
setwd("CHOL")

# Download files
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S2211124717302140-mmc1.pdf"
download.file(curl, "document_s1.pdf")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S2211124717302140-mmc2.xlsx"
download.file(curl, "tables.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S2211124717302140-mmc3.pdf"
download.file(curl, "ResourceIntegrative Genomic Analysis of Cholangiocarcinoma Identifies Distinct IDH-Mutant Molecular Profiles.pdf")

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

setwd("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/reference_data")
file_exists <- file.exists("COAD")
if (file_exists == FALSE) {
  dir.create("COAD")
}
setwd("COAD")

# Download files
curl <- "https://static-content.springer.com/esm/art%3A10.1038%2Fnature11252/MediaObjects/41586_2012_BFnature11252_MOESM70_ESM.zip"
download.file(curl, "supplementary_table2.zip")
curl <- "https://static-content.springer.com/esm/art%3A10.1038%2Fnature11252/MediaObjects/41586_2012_BFnature11252_MOESM71_ESM.zip"
download.file(curl, "supplementary_table1_and_3_12.zip")
curl <- "https://www.nature.com/articles/nature11252.pdf"
download.file(curl, "Comprehensive molecular characterization of human colon and rectal cancer.pdf")

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

setwd("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/reference_data")
file_exists <- file.exists("ESCA")
if (file_exists == FALSE) {
  dir.create("ESCA")
}
setwd("ESCA")

# Download files
curl <- "https://static-content.springer.com/esm/art%3A10.1038%2Fnature20805/MediaObjects/41586_2017_BFnature20805_MOESM91_ESM.xlsx"
download.file(curl, "table_s1.xlsx")
curl <- "https://static-content.springer.com/esm/art%3A10.1038%2Fnature20805/MediaObjects/41586_2017_BFnature20805_MOESM92_ESM.xlsx"
download.file(curl, "table_s2.xlsx")
curl <- "https://static-content.springer.com/esm/art%3A10.1038%2Fnature20805/MediaObjects/41586_2017_BFnature20805_MOESM93_ESM.xlsx"
download.file(curl, "table_s3.xlsx")
curl <- "https://static-content.springer.com/esm/art%3A10.1038%2Fnature20805/MediaObjects/41586_2017_BFnature20805_MOESM94_ESM.xlsx"
download.file(curl, "table_s4.xlsx")
curl <- "https://static-content.springer.com/esm/art%3A10.1038%2Fnature20805/MediaObjects/41586_2017_BFnature20805_MOESM95_ESM.xlsx"
download.file(curl, "table_s5.xlsx")
curl <- "https://static-content.springer.com/esm/art%3A10.1038%2Fnature20805/MediaObjects/41586_2017_BFnature20805_MOESM96_ESM.xlsx"
download.file(curl, "table_s6.xlsx")
curl <- "https://static-content.springer.com/esm/art%3A10.1038%2Fnature20805/MediaObjects/41586_2017_BFnature20805_MOESM97_ESM.xlsx"
download.file(curl, "table_s7.xlsx")
curl <- "https://static-content.springer.com/esm/art%3A10.1038%2Fnature20805/MediaObjects/41586_2017_BFnature20805_MOESM98_ESM.pdf"
download.file(curl, "supplementary_information.pdf")
curl <- "https://www.nature.com/articles/nature20805.pdf"
download.file(curl, "Integrated genomic characterization of oesophageal carcinoma.pdf")

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

setwd("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/reference_data")
file_exists <- file.exists("GBM")
if (file_exists == FALSE) {
  dir.create("GBM")
}
setwd("GBM")

# Download files
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867413012087-mmc1.xls"
download.file(curl, "table_s1.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867413012087-mmc2.xlsx"
download.file(curl, "table_s2.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867413012087-mmc3.xlsx"
download.file(curl, "table_s3.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867413012087-mmc4.xlsx"
download.file(curl, "table_s4.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867413012087-mmc5.xlsx"
download.file(curl, "table_s5.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867413012087-mmc6.xlsx"
download.file(curl, "table_s6.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867413012087-mmc7.xlsx"
download.file(curl, "table_s7.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867413012087-mmc8.pdf"
download.file(curl, "document_s1.pdf")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867413012087-mmc9.pdf"
download.file(curl, "The Somatic Genomic Landscapeof Glioblastoma.pdf")

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

setwd("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/reference_data")
file_exists <- file.exists("HNSC")
if (file_exists == FALSE) {
  dir.create("HNSC")
}
setwd("HNSC")

# Download files
curl <- "https://static-content.springer.com/esm/art%3A10.1038%2Fnature14129/MediaObjects/41586_2015_BFnature14129_MOESM116_ESM.zip"
download.file(curl, "supplementary_data.zip")
curl <- "https://www.nature.com/articles/nature14129.pdf"
download.file(curl, "Comprehensive genomic characterization of head and neck squamous cell carcinomas.pdf")

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

setwd("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/reference_data")
file_exists <- file.exists("KICH")
if (file_exists == FALSE) {
  dir.create("KICH")
}
setwd("KICH")

# Download files
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610814003043-mmc1.pdf"
download.file(curl, "document_s1.pdf")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610814003043-mmc2.xlsx"
download.file(curl, "table_s1.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610814003043-mmc3.xlsx"
download.file(curl, "table_s2.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610814003043-mmc4.xlsx"
download.file(curl, "table_s3.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610814003043-mmc5.xlsx"
download.file(curl, "table_s4.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610814003043-mmc6.xlsx"
download.file(curl, "table_s5.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610814003043-mmc7.xlsx"
download.file(curl, "table_s6.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610814003043-mmc8.xlsx"
download.file(curl, "table_s7.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610814003043-mmc9.xls"
download.file(curl, "table_s8.xls")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610814003043-mmc10.pdf"
download.file(curl, "The Somatic Genomic Landscape of Chromophobe Renal Cell Carcinoma.pdf")

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

setwd("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/reference_data")
file_exists <- file.exists("KIRC")
if (file_exists == FALSE) {
  dir.create("KIRC")
}
setwd("KIRC")

# Download files
curl <- "https://static-content.springer.com/esm/art%3A10.1038%2Fnature12222/MediaObjects/41586_2013_BFnature12222_MOESM35_ESM.zip"
download.file(curl, "supplementary_data.zip")
curl <- "https://www.nature.com/articles/nature12222.pdf"
download.file(curl, "Comprehensive molecular characterization of clear cell renal cell carcinoma.pdf")

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

setwd("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/reference_data")
file_exists <- file.exists("KIRP")
if (file_exists == FALSE) {
  dir.create("KIRP")
}
setwd("KIRP")

# Download files
curl <- "https://www.nejm.org/doi/suppl/10.1056/NEJMoa1505917/suppl_file/nejmoa1505917_appendix_1.pdf"
download.file(curl, "supplementary_appendix1.pdf")
curl <- "https://www.nejm.org/doi/suppl/10.1056/NEJMoa1505917/suppl_file/nejmoa1505917_appendix_2.xlsx"
download.file(curl, "supplementary_appendix2.xlsx")
curl <- "https://www.nejm.org/doi/suppl/10.1056/NEJMoa1505917/suppl_file/nejmoa1505917_appendix_3.xlsx"
download.file(curl, "supplementary_appendix3.xlsx")
curl <- "https://www.nejm.org/doi/suppl/10.1056/NEJMoa1505917/suppl_file/nejmoa1505917_appendix_4.xlsx"
download.file(curl, "supplementary_appendix4.xlsx")
curl <- "https://www.nejm.org/doi/suppl/10.1056/NEJMoa1505917/suppl_file/nejmoa1505917_appendix_5.xlsx"
download.file(curl, "supplementary_appendix5.xlsx")
curl <- "https://www.nejm.org/doi/suppl/10.1056/NEJMoa1505917/suppl_file/nejmoa1505917_appendix_6.xlsx"
download.file(curl, "supplementary_appendix6.xlsx")
curl <- "https://www.nejm.org/doi/suppl/10.1056/NEJMoa1505917/suppl_file/nejmoa1505917_appendix_7.xlsx"
download.file(curl, "supplementary_appendix7.xlsx")
curl <- "https://www.nejm.org/doi/suppl/10.1056/NEJMoa1505917/suppl_file/nejmoa1505917_appendix_8.xlsx"
download.file(curl, "supplementary_appendix8.xlsx")
curl <- "https://www.nejm.org/doi/suppl/10.1056/NEJMoa1505917/suppl_file/nejmoa1505917_appendix_9.xlsx"
download.file(curl, "supplementary_appendix9.xlsx")
curl <- "https://www.nejm.org/doi/suppl/10.1056/NEJMoa1505917/suppl_file/nejmoa1505917_appendix_10.xlsx"
download.file(curl, "supplementary_appendix10.xlsx")
curl <- "https://www.nejm.org/doi/pdf/10.1056/NEJMoa1505917?articleTools=true"
download.file(curl, "Comprehensive Molecular Characterization of Papillary Renal-Cell Carcinoma.pdf")

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

setwd("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/reference_data")
file_exists <- file.exists("LAML")
if (file_exists == FALSE) {
  dir.create("LAML")
}
setwd("LAML")

# Download files
curl <- "https://cbioportal-datahub.s3.amazonaws.com/laml_tcga_pub.tar.gz"
download.file(curl, "aml_nejm_2013.tar")
curl <- "https://www.nejm.org/doi/pdf/10.1056/NEJMoa1301689?articleTools=true"
download.file(curl, "Genomic and epigenomic landscapes of adult de novo acute myeloid leukemia.pdf")

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

setwd("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/reference_data")
file_exists <- file.exists("LGG")
if (file_exists == FALSE) {
  dir.create("LGG")
}
setwd("LGG")

# Download files
curl <- "https://www.nejm.org/doi/suppl/10.1056/NEJMoa1402121/suppl_file/nejmoa1402121_appendix_1.pdf"
download.file(curl, "supplementary_appendix1.pdf")
curl <- "https://www.nejm.org/doi/suppl/10.1056/NEJMoa1402121/suppl_file/nejmoa1402121_appendix_2.xlsx"
download.file(curl, "supplementary_appendix2.xlsx")
curl <- "https://www.nejm.org/doi/suppl/10.1056/NEJMoa1402121/suppl_file/nejmoa1402121_appendix_3.xlsx"
download.file(curl, "supplementary_appendix3.xlsx")
curl <- "https://www.nejm.org/doi/suppl/10.1056/NEJMoa1402121/suppl_file/nejmoa1402121_appendix_4.xlsx"
download.file(curl, "supplementary_appendix4.xlsx")
curl <- "https://www.nejm.org/doi/suppl/10.1056/NEJMoa1402121/suppl_file/nejmoa1402121_appendix_5.xlsx"
download.file(curl, "supplementary_appendix5.xlsx")
curl <- "https://www.nejm.org/doi/suppl/10.1056/NEJMoa1402121/suppl_file/nejmoa1402121_appendix_6.xlsx"
download.file(curl, "supplementary_appendix6.xlsx")
curl <- "https://www.nejm.org/doi/suppl/10.1056/NEJMoa1402121/suppl_file/nejmoa1402121_appendix_7.xlsx"
download.file(curl, "supplementary_appendix7.xlsx")
curl <- "https://www.nejm.org/doi/suppl/10.1056/NEJMoa1402121/suppl_file/nejmoa1402121_appendix_8.xls"
download.file(curl, "supplementary_appendix8.xls")
curl <- "https://www.nejm.org/doi/pdf/10.1056/NEJMoa1402121?articleTools=true"
download.file(curl, "Comprehensive, Integrative Genomic Analysis of Diffuse Lower-Grade Gliomas.pdf")

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

setwd("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/reference_data")
file_exists <- file.exists("LIHC")
if (file_exists == FALSE) {
  dir.create("LIHC")
}
setwd("LIHC")

# Download files
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867417306396-mmc1.xlsx"
download.file(curl, "table_s1.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867417306396-mmc2.xlsx"
download.file(curl, "table_s2.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867417306396-mmc3.xlsx"
download.file(curl, "table_s3.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867417306396-mmc4.xlsx"
download.file(curl, "table_s4.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867417306396-mmc5.xlsx"
download.file(curl, "table_s5.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867417306396-mmc6.xlsx"
download.file(curl, "table_s6.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867417306396-mmc7.xlsx"
download.file(curl, "table_s7.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867417306396-mmc8.xlsx"
download.file(curl, "table_s8.xlsx")
curl <- "https://www.sciencedirect.com/science/article/pii/S0092867417306396/pdfft?isDTMRedir=true&download=true"
download.file(curl, "Comprehensive and Integrative Genomic Characterization of Hepatocellular Carcinoma.pdf")

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

setwd("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/reference_data")
file_exists <- file.exists("LUAD")
if (file_exists == FALSE) {
  dir.create("LUAD")
}
setwd("LUAD")

# Download files
curl <- "https://static-content.springer.com/esm/art%3A10.1038%2Fnature13385/MediaObjects/41586_2014_BFnature13385_MOESM21_ESM.xlsx"
download.file(curl, "supplementary_tables.xlsx")
curl <- "https://www.nature.com/articles/nature13385.pdf"
download.file(curl, "Comprehensive molecular profiling of lung adenocarcinoma.pdf")

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

setwd("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/reference_data")
file_exists <- file.exists("LUSC")
if (file_exists == FALSE) {
  dir.create("LUSC")
}
setwd("LUSC")

# Download files
curl <- "https://static-content.springer.com/esm/art%3A10.1038%2Fnature11404/MediaObjects/41586_2012_BFnature11404_MOESM286_ESM.zip"
download.file(curl, "supplementary_data.zip")
curl <- "https://www.nature.com/articles/nature11404.pdf"
download.file(curl, "Comprehensive genomic characterization of squamous cell lung cancers.pdf")

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

setwd("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/reference_data")
file_exists <- file.exists("MESO")
if (file_exists == FALSE) {
  dir.create("MESO")
}
setwd("MESO")

# Download files
curl <- "https://cancerdiscovery.aacrjournals.org/highwire/filestream/44394/field_highwire_adjunct_files/0/205173_2_supp_5073238_pg5mpy.pdf"
download.file(curl, "supplementary_figures.pdf")
curl <- "https://cancerdiscovery.aacrjournals.org/highwire/filestream/44394/field_highwire_adjunct_files/1/205173_2_supp_5073239_pg14sl.xlsx"
download.file(curl, "table_s1.xlsx")
curl <- "https://cancerdiscovery.aacrjournals.org/highwire/filestream/44394/field_highwire_adjunct_files/2/205173_2_supp_5073240_pg14sl.xlsx"
download.file(curl, "table_s2.xlsx")
curl <- "https://cancerdiscovery.aacrjournals.org/highwire/filestream/44394/field_highwire_adjunct_files/3/205173_2_supp_5073241_pg14sl.docx"
download.file(curl, "table_s3.docx")
curl <- "https://cancerdiscovery.aacrjournals.org/highwire/filestream/44394/field_highwire_adjunct_files/4/205173_2_supp_5073242_pg14sl.xlsx"
download.file(curl, "table_s4.xlsx")
curl <- "https://cancerdiscovery.aacrjournals.org/highwire/filestream/44394/field_highwire_adjunct_files/5/205173_2_supp_5073243_pggmgr.xlsx"
download.file(curl, "table_s5.xlsx")
curl <- "https://cancerdiscovery.aacrjournals.org/highwire/filestream/44394/field_highwire_adjunct_files/6/205173_2_supp_5073244_pg14sl.xlsx"
download.file(curl, "table_s6.xlsx")
curl <- "https://cancerdiscovery.aacrjournals.org/highwire/filestream/44394/field_highwire_adjunct_files/7/205173_2_supp_5073245_pg14sl.xlsx"
download.file(curl, "table_s7.xlsx")
curl <- "https://cancerdiscovery.aacrjournals.org/content/8/12/1548.full-text.pdf"
download.file(curl, "Integrative Molecular Characterization of Malignant Pleural Mesothelioma.pdf")

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

setwd("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/reference_data")
file_exists <- file.exists("OV")
if (file_exists == FALSE) {
  dir.create("OV")
}
setwd("OV")

# Download files
curl <- "https://static-content.springer.com/esm/art%3A10.1038%2Fnature10166/MediaObjects/41586_2011_BFnature10166_MOESM44_ESM.zip"
download.file(curl, "supplementary_data.zip")
curl <- "https://www.nature.com/articles/nature10166.pdf"
download.file(curl, "Integrated genomic analyses of ovarian carcinoma.pdf")

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

setwd("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/reference_data")
file_exists <- file.exists("PAAD")
if (file_exists == FALSE) {
  dir.create("PAAD")
}
setwd("PAAD")

# Download files
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610817302994-mmc1.pdf"
download.file(curl, "document_s1.pdf")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610817302994-mmc2.xlsx"
download.file(curl, "table_s1.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610817302994-mmc3.xlsx"
download.file(curl, "table_s2.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610817302994-mmc4.xlsx"
download.file(curl, "table_s3.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610817302994-mmc5.xlsx"
download.file(curl, "table_s4.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610817302994-mmc6.xlsx"
download.file(curl, "table_s5.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610817302994-mmc7.xlsx"
download.file(curl, "table_s6.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610817302994-mmc8.xlsx"
download.file(curl, "table_s7.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610817302994-mmc9.xlsx"
download.file(curl, "table_s8.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610817302994-mmc10.xlsx"
download.file(curl, "table_s9.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610817302994-mmc11.pdf"
download.file(curl, "Integrated Genomic Characterization of PancreaticDuctal Adenocarcinoma.pdf")

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

setwd("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/reference_data")
file_exists <- file.exists("PCPG")
if (file_exists == FALSE) {
  dir.create("PCPG")
}
setwd("PCPG")

# Download files
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610817300016-mmc1.pdf"
download.file(curl, "document_s1.pdf")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610817300016-mmc2.xlsx"
download.file(curl, "table_s1.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610817300016-mmc3.xls"
download.file(curl, "table_s2.xls")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610817300016-mmc4.pdf"
download.file(curl, "Comprehensive Molecular Characterization ofPheochromocytoma and Paraganglioma.pdf")

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

setwd("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/reference_data")
file_exists <- file.exists("PRAD")
if (file_exists == FALSE) {
  dir.create("PRAD")
}
setwd("PRAD")

# Download files
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867415013392-mmc1.pdf"
download.file(curl, "document_s1.pdf")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867415013392-mmc2.xls"
download.file(curl, "table_s1.xls")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867415013392-mmc3.xls"
download.file(curl, "table_s2.xls")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867415013392-mmc4.xls"
download.file(curl, "table_s3.xls")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867415013392-mmc5.pdf"
download.file(curl, "The Molecular Taxonomy of Primary Prostate Cancer.pdf")

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

setwd("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/reference_data")
file_exists <- file.exists("SARC_Adult")
if (file_exists == FALSE) {
  dir.create("SARC_Adult")
}
setwd("SARC_Adult")

# Download files
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867417312035-mmc1.xlsx"
download.file(curl, "table_s1.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867417312035-mmc2.pdf"
download.file(curl, "table_s2.pdf")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867417312035-mmc3.xlsx"
download.file(curl, "table_s3.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867417312035-mmc4.xlsx"
download.file(curl, "table_s4.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867417312035-mmc5.xlsx"
download.file(curl, "table_s5.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867417312035-mmc6.xlsx"
download.file(curl, "table_s6.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867417312035-mmc7.xls"
download.file(curl, "table_s7.xls")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867417312035-mmc8.pdf"
download.file(curl, "methods_s1.pdf")
curl <- "https://www.sciencedirect.com/science/article/pii/S0092867417312035/pdfft?isDTMRedir=true&download=true"
download.file(curl, "Comprehensive and Integrated GenomicCharacterization of Adult Soft Tissue Sarcomas.pdf")

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

setwd("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/reference_data")
file_exists <- file.exists("SKCM")
if (file_exists == FALSE) {
  dir.create("SKCM")
}
setwd("SKCM")

# Download files
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867415006340-mmc1.pdf"
download.file(curl, "document_s1.pdf")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867415006340-mmc2.xlsx"
download.file(curl, "table_s1.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867415006340-mmc3.xlsx"
download.file(curl, "table_s2.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867415006340-mmc4.xlsx"
download.file(curl, "table_s3.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867415006340-mmc5.xlsx"
download.file(curl, "table_s4.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867415006340-mmc6.pdf"
download.file(curl, "Genomic Classification of Cutaneous Melanoma.pdf")

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

setwd("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/reference_data")
file_exists <- file.exists("STAD")
if (file_exists == FALSE) {
  dir.create("STAD")
}
setwd("STAD")

# Download files
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610818301144-mmc1.pdf"
download.file(curl, "document_s1.pdf")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610818301144-mmc2.xlsx"
download.file(curl, "table_s1.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610818301144-mmc3.xlsx"
download.file(curl, "table_s2.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610818301144-mmc4.xlsx"
download.file(curl, "table_s3.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610818301144-mmc5.xlsx"
download.file(curl, "table_s4.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610818301144-mmc6.xlsx"
download.file(curl, "table_s5.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610818301144-mmc7.xlsx"
download.file(curl, "table_s6.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610818301144-mmc8.xlsx"
download.file(curl, "table_s7.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610818301144-mmc9.pdf"
download.file(curl, "Comparative Molecular Analysis of Gastrointestinal Adenocarcinomas.pdf")

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

setwd("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/reference_data")
file_exists <- file.exists("TGCT")
if (file_exists == FALSE) {
  dir.create("TGCT")
}
setwd("TGCT")

# Download files
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S221112471830785X-mmc1.pdf"
download.file(curl, "document_s1_table_s2.pdf")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S221112471830785X-mmc2.xlsx"
download.file(curl, "table_s1.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S221112471830785X-mmc3.xlsx"
download.file(curl, "table_s3.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S221112471830785X-mmc4.xlsx"
download.file(curl, "table_s4.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S221112471830785X-mmc5.xlsx"
download.file(curl, "table_s5.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S221112471830785X-mmc6.xlsx"
download.file(curl, "table_s6.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S221112471830785X-mmc7.xlsx"
download.file(curl, "table_s7.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S221112471830785X-mmc8.pdf"
download.file(curl, "Integrated Molecular Characterization of TesticularGerm Cell Tumors.pdf")

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

setwd("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/reference_data")
file_exists <- file.exists("THCA_PTC")
if (file_exists == FALSE) {
  dir.create("THCA_PTC")
}
setwd("THCA_PTC")

# Download files
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867414012380-mmc1.pdf"
download.file(curl, "document_s1.pdf")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867414012380-mmc2.pdf"
download.file(curl, "document_s2.pdf")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867414012380-mmc3.xlsx"
download.file(curl, "table_s2.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867414012380-mmc4.xlsx"
download.file(curl, "table_s4.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867414012380-mmc5.xlsx"
download.file(curl, "table_s6.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867414012380-mmc6.zip"
download.file(curl, "data_s1.zip")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867414012380-mmc7.pdf"
download.file(curl, "document_s3.pdf")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867414012380-mmc8.pdf"
download.file(curl, "Integrated Genomic Characterizationof Papillary Thyroid Carcinoma.pdf")

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

setwd("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/reference_data")
file_exists <- file.exists("THYM_TETs")
if (file_exists == FALSE) {
  dir.create("THYM_TETs")
}
setwd("THYM_TETs")

# Download files
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610818300035-mmc1.pdf"
download.file(curl, "document_s1.pdf")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610818300035-mmc2.xlsx"
download.file(curl, "table_s1.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610818300035-mmc3.pdf"
download.file(curl, "The Integrated Genomic Landscape of ThymicEpithelial Tumors.pdf")

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

setwd("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/reference_data")
file_exists <- file.exists("UCEC")
if (file_exists == FALSE) {
  dir.create("UCEC")
}
setwd("UCEC")

# Download files
curl <- "https://static-content.springer.com/esm/art%3A10.1038%2Fnature12113/MediaObjects/41586_2013_BFnature12113_MOESM89_ESM.zip"
download.file(curl, "supplementary_data.zip")
curl <- "https://www.nature.com/articles/nature12113.pdf"
download.file(curl, "Integrated genomic characterization of endometrial carcinoma.pdf")

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

setwd("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/reference_data")
file_exists <- file.exists("UCS")
if (file_exists == FALSE) {
  dir.create("UCS")
}
setwd("UCS")

# Download files
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610817300533-mmc1.pdf"
download.file(curl, "document_s1.pdf")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610817300533-mmc2.xlsx"
download.file(curl, "table_s1.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610817300533-mmc3.xlsx"
download.file(curl, "table_s2.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610817300533-mmc4.xlsx"
download.file(curl, "table_s3.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610817300533-mmc5.xlsx"
download.file(curl, "table_s4.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610817300533-mmc6.xlsx"
download.file(curl, "table_s5.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610817300533-mmc7.xlsx"
download.file(curl, "table_s6.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610817300533-mmc8.xlsx"
download.file(curl, "table_s7.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610817300533-mmc9.pdf"
download.file(curl, "Integrated Molecular Characterization of UterineCarcinosarcoma.pdf")

# Environment Clear
rm(list = ls(all.names = TRUE))
ls(all.names = TRUE)

setwd("/Volumes/G_DRIVEmobile/Senn_CDK_Cyclins_Cancer/reference_data")
file_exists <- file.exists("UVM")
if (file_exists == FALSE) {
  dir.create("UVM")
}
setwd("UVM")

# Download files
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610817302957-mmc1.pdf"
download.file(curl, "document_s1.pdf")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610817302957-mmc2.xlsx"
download.file(curl, "table_s1.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610817302957-mmc3.xlsx"
download.file(curl, "table_s2.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610817302957-mmc4.xlsx"
download.file(curl, "table_s3.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610817302957-mmc5.xlsx"
download.file(curl, "table_s4.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610817302957-mmc6.xlsx"
download.file(curl, "table_s5.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610817302957-mmc7.xlsx"
download.file(curl, "table_s6.xlsx")
curl <- "https://ars.els-cdn.com/content/image/1-s2.0-S1535610817302957-mmc8.pdf"
download.file(curl, "ArticleIntegrative Analysis Identifies Four Molecular andClinical Subsets in Uveal Melanoma.pdf")

traceback()

