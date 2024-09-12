library(tidyverse)
library(stringi)
library(FactoMineR)
library(factoextra)

# Read vcf files -------------------------------------------------------------
extract_AIMs <- function(directory, AIMs, samples) {
  
  AIMs_rsids <- AIMs
  filepaths <- file.path(directory, samples$filename)
  file_names <- basename(filepaths)
  file_names_alias <- paste0("File", c(1:length(file_names)))
  names(file_names_alias) <- file_names
  df <- matrix(nrow = length(filepaths), ncol = length(AIMs_rsids)) %>% 
    as.data.frame(row.names = file_names_alias)
  colnames(df) <- AIMs_rsids
  
  print(paste("Total files", "=", length(file_names_alias)))
  
  for (filepath in filepaths) {
    
    vcf_file <- read.csv(file = filepath, 
                     header = FALSE, 
                     comment.char = "#", 
                     sep = "\t")
    colnames(vcf_file) <- c("rsid", "chromosome", "position", "genotype")
    vcf_file <- vcf_file[!(vcf_file$rsid %>% duplicated()), ]
    rownames(vcf_file) <- vcf_file$rsid
    genotypes <- vcf_file[AIMs_rsids, "genotype"]
    genotypes[is.na(genotypes)] <- "--"
    df[file_names_alias[basename(filepath)], ] <- genotypes
    print(file_names_alias[basename(filepath)])
    
  }
  
  return(list(df, file_names_alias))
  
}

# Count nucleotides ------------------------------------------------------------
genotype_to_numeric <- function(sample_data) {
  
  SNP_names <- colnames(sample_data)
  SNP_names_alias <- paste0("SNP", c(1:length(SNP_names)))
  names(SNP_names_alias) <- SNP_names
  df <- matrix(nrow = dim(sample_data)[1], 
               ncol = dim(sample_data)[2]*4) %>% 
    as.data.frame(row.names = rownames(sample_data))
  col_names <- vector(mode = "character", length = dim(sample_data)[2]*4)
  i <- 1
  for (SNP_name_alias in SNP_names_alias) {
    
    col_names[c(i : (i + 4 - 1))] <- paste0(SNP_name_alias, c("A", "T", "G", "C"))
    i <- i + 4
  }
  colnames(df) <- col_names
  for (file in rownames(df)) {
    
    signature <- sample_data[file, ]
    
    counts_nuc <- stri_count(str = signature, fixed = "A")
    df[file, seq(1, 4*dim(sample_data)[2], 4)] <- counts_nuc
    counts_nuc <- stri_count(str = signature, fixed = "T")
    df[file, seq(2, 4*dim(sample_data)[2], 4)] <- counts_nuc
    counts_nuc <- stri_count(str = signature, fixed = "G")
    df[file, seq(3, 4*dim(sample_data)[2], 4)] <- counts_nuc
    counts_nuc <- stri_count(str = signature, fixed = "C")
    df[file, seq(4, 4*dim(sample_data)[2], 4)] <- counts_nuc
    
  }
  
  return(list(df, SNP_names_alias))
  
}

# Analysis of test folder ------------------------------------------------------
AIMs <- read_lines(file = "/media/singhlab/FF39-032F/Final_project/GWAS_script/3_Association_GWAS/assoc_results_adjusted_trimmed.txt", skip = 1)
AIMs <- sub(".*\\brs([0-9]+)\\b.*", "rs\\1", AIMs)
input_samples_GWAS <- read.table(file = "/media/singhlab/FF39-032F/Final_project/data_integrated/samples_used_for_GWAS.txt", sep = ";", header = TRUE)
output_samples_GWAS <- read.table(file = "/media/singhlab/FF39-032F/Final_project/GWAS_script/3_Association_GWAS/23andme_TR_13.fam", sep = " ")
samples_GWAS <- input_samples_GWAS[input_samples_GWAS$FID %in% output_samples_GWAS$V1, ]
directory <- "/media/singhlab/FF39-032F/Final_project/data_raw"
extract_AIMs_out <- extract_AIMs(directory = directory, AIMs = AIMs, samples = samples_GWAS)
sample_data <- extract_AIMs_out[[1]]
file_names_alias <- extract_AIMs_out[[2]]
genotype_to_numeric_out <- genotype_to_numeric(sample_data)
sample_data_numeric <- genotype_to_numeric_out[[1]]
SNP_names_alias <- genotype_to_numeric_out[[2]]

samples_GWAS_phenotypes <- samples_GWAS$Binary_category %>% as.data.frame()
colnames(samples_GWAS_phenotypes) <- "Phenotypes"
sample_data_numeric <- cbind(samples_GWAS_phenotypes, sample_data_numeric)

write.table(x = sample_data_numeric, file = "/media/singhlab/FF39-032F/Final_project/out_ML/data.csv", quote = FALSE, sep = ",")
SNP_names <- data.frame(Alias = SNP_names_alias, row.names = names(SNP_names_alias))
file_names <- data.frame(Alisa = file_names_alias, row.names = names(file_names_alias))
write.table(SNP_names, file = "/media/singhlab/FF39-032F/Final_project/out_ML/SNP_names_dict.csv", quote = FALSE, sep = ",")
write.table(file_names, file = "/media/singhlab/FF39-032F/Final_project/out_ML/file_names_dict.csv", quote = FALSE, sep = ",")
