library(tidyverse)
library(rstudioapi)
library(stringi)

df_transpose <- function(df){
  col_names <- rownames(df)
  row_names <- colnames(df)
  df <- as.matrix(df)
  df <- t(df)
  df <- as.data.frame(df, row.names = row_names)
  colnames(df) <- col_names
  return(df)
}

current_dir <- getActiveDocumentContext()$path %>% dirname()
samples_data <- read.csv(file = file.path(current_dir, "phenotypes_202304020100.csv"), 
                         header = TRUE, sep = ";")
samples_data_used <- samples_data[, c("user_id", "genotype_filename", "date_of_birth", "chrom_sex", "Tongue.roller")]
samples_data_used <- samples_data_used[!(samples_data_used$Tongue.roller == "-"), ]
categories <- unique(samples_data_used$Tongue.roller)
print(categories)
dict <- vector(mode = "integer", length = length(categories))
dict <- c(1, 2, 1, 2, 2, 1, 2, 1, 1, 1, 1, 2)
names(dict) <- categories
samples_data_used$Binary_category <- dict[samples_data_used$Tongue.roller]
samples_data_used$date_of_birth[samples_data_used$date_of_birth == "rather not say"] <- "unknown"
samples_data_used$chrom_sex[samples_data_used$chrom_sex == "rather not say"] <- "unknown"
temp2 <- stri_split_fixed(str = samples_data_used$genotype_filename, 
                          pattern = ".") %>% 
  as.data.frame() %>% 
  df_transpose()
colnames(temp2) <- c("user_id", "filetype", "file_number")
samples_data_used <- cbind(samples_data_used, temp2)
samples_data_used$filename <- paste0("user", 
                                     samples_data_used$user_id, 
                                     "_file", 
                                     samples_data_used$file_number, 
                                     "_yearofbirth_", 
                                     samples_data_used$date_of_birth, 
                                     "_sex_", 
                                     samples_data_used$chrom_sex, 
                                     ".", 
                                     samples_data_used$filetype, 
                                     ".txt")

samples_data_used <- samples_data_used[(samples_data_used$filetype == "23andme"), ]
samples_data_used$size <- file.info(file.path(current_dir, "data_raw", samples_data_used$filename))$size

setwd(current_dir)
opensnp_urls <- file.path("https://opencravat.org/data/open-snp", 
                          samples_data_used$filename)
opensnp_dests <- file.path(current_dir, "data_raw", samples_data_used$filename)

for (i in c(1:length(samples_data_used$filename))) {
  
  download.file(url = opensnp_urls[i], destfile = opensnp_dests[i], method = "wget")
  
}

# Remove duplicates ------------------------------------------------------------
duplicated_user_ids <- samples_data_used[duplicated(samples_data_used$user_id), ]$user_id %>% unique()
for (user_id in duplicated_user_ids) {
  
  rep_files <- samples_data_used[samples_data_used$user_id == user_id, ]
  max_loc <- rep_files$size %>% which.max()
  max_file_name <- rep_files$filename[max_loc]
  remove_file_names <- rep_files$filename[!(rep_files$filename %in% max_file_name)]
  samples_data_used <- samples_data_used[!(samples_data_used$filename %in% remove_file_names), ]
  
}

# Catch unreadable files -------------------------------------------------------
files_list <- file.path(current_dir, "data_raw", samples_data_used$filename)
error_files <- list()
for (file_name in files_list) {
  
  read_state <- tryCatch(
    {
      read_lines(file = file_name, n_max = 1)
      #read_tsv(file = file_name, comment = "#", col_names = FALSE)
    }, 
    error = function(e) e
  )
  if (inherits(read_state, "error")) {
    error_files[[length(error_files) + 1]] <- file_name
  }
  
}
error_files <- unlist(error_files)
samples_data_used <- samples_data_used[!(samples_data_used$filename %in% basename(error_files)), ]

# Segregate builds -------------------------------------------------------------
# some 23andme files are actually VCF. These files are also unintentionally 
# filtered in this step
files_list <- file.path(current_dir, "data_raw", samples_data_used$filename)
build_info <- vector(mode = "character", length = length(files_list))
names(build_info) <- files_list
for (file_name in files_list) {
  
  ref_genome <- read_lines(file = file_name, skip = 11, n_max = 1)
  ref_genome <- gsub(".*build\\s+(\\d+).*", "\\1", ref_genome)
  build_info[file_name] <- ref_genome
  
}
samples_data_used <- samples_data_used[build_info == "37", ]

# Create binary dataset --------------------------------------------------------
individual_files <- samples_data_used$filename
individual_phenotypes <- samples_data_used$Binary_category
for (n_file in c(1:length(individual_files))) {
  
  in_file <- file.path(current_dir, "data_raw", individual_files[n_file])
  file_name_without_ext <- sub(pattern = "(.*?)\\..*$", replacement = "\\1", individual_files[n_file])
  out_file <- file.path(current_dir, "data_binary", file_name_without_ext)
  system(command = paste("~/plink-1.07-x86_64/plink", 
                         "--23file", 
                         in_file, 
                         n_file, 
                         n_file, 
                         "0", 
                         individual_phenotypes[n_file], 
                         "0", 
                         "0", 
                         "--out", 
                         out_file,  
                         sep = " "))
  print(n_file)
   
}

# Remove erroneous files -------------------------------------------------------
log_files <- list.files(path = file.path(current_dir, "data_binary"), pattern = ".*.log$", full.names = TRUE)
error_files <- list()
for (log_file in log_files) {
  
  file_data <- readLines(log_file)
  error_line_loc <- grep(pattern = "^Error:", file_data)
  if (!(length(error_line_loc) == 0)) {
    file_loc_txt <- file_data[error_line_loc + 1]
    error_msg <- file_data[error_line_loc + 2]
    error_files[[basename(log_file)]] <- c(basename(log_file), file_loc_txt, error_msg)
  }
  
}

# Correcting Chromosome order for File user10632_file8882_yearofbirth_1958_sex_XY.23andme.txt
temp <- read.table(file = "/media/singhlab/FF39-032F/Final_project/data_raw/user10632_file8882_yearofbirth_1958_sex_XY.23andme.txt", comment.char = "#", header = FALSE)
unique(temp$V2)
# Output: [1] "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "MT" "X"  "Y" 
# Need to convert it to: [1] "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X"  "Y"  "MT"
temp_MT <- temp[temp$V2 == "MT", ]
temp <- temp[!(temp$V2 == "MT"), ]
temp <- rbind(temp, temp_MT)
unique(temp$V2)
# Output: [1] "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X"  "Y"  "MT"
write.table(x = temp, file = "/media/singhlab/FF39-032F/Final_project/data_raw/user10632_file8882_yearofbirth_1958_sex_XY.23andme_temp.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Correcting Chromosome order for File user4577_file8223_yearofbirth_unknown_sex_unknown.23andme.txt
temp <- read.table(file = "/media/singhlab/FF39-032F/Final_project/data_raw/user4577_file8223_yearofbirth_unknown_sex_unknown.23andme.txt", comment.char = "#", header = FALSE)
unique(temp$V2)
# Output: [1] "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "MT" "X"  "Y" 
# Need to convert it to: [1] "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X"  "Y"  "MT"
temp_MT <- temp[temp$V2 == "MT", ]
temp <- temp[!(temp$V2 == "MT"), ]
temp <- rbind(temp, temp_MT)
unique(temp$V2)
# Output: [1] "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X"  "Y"  "MT"
write.table(x = temp, file = "/media/singhlab/FF39-032F/Final_project/data_raw/user4577_file8223_yearofbirth_unknown_sex_unknown.23andme_temp.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

error_files <- error_files %>% as.data.frame() %>% df_transpose()
error_files <- sub(pattern = "(.*?)\\..*$", replacement = "\\1", rownames(error_files))
error_files <- paste(error_files, "23andme", "txt", sep = ".")
samples_data_used <- samples_data_used[!(samples_data_used$filename %in% error_files), ]

# Create integrated binary dataset --------------------------------------------- 
complete_valid_binary_files <- sub(pattern = "(.*?)\\..*$", replacement = "\\1", samples_data_used$filename)
bed_files <- file.path(current_dir, "data_binary", paste0(complete_valid_binary_files, ".bed"))
bim_files <- file.path(current_dir, "data_binary", paste0(complete_valid_binary_files, ".bim"))
fam_files <- file.path(current_dir, "data_binary", paste0(complete_valid_binary_files, ".fam"))
df_filenames <- data.frame(bed = bed_files, bim = bim_files, fam = fam_files)
first_file <- file.path(current_dir, "data_binary", complete_valid_binary_files[1])
df_filenames <- df_filenames[c(-1), ]
write.table(x = df_filenames, file = file.path(current_dir, "data_integrated", "merge_files.txt"), sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)
written_file <- file.path(current_dir, "data_integrated", "merge_files.txt")
out_file <- file.path(current_dir, "data_integrated", "23andme_TR")
system(command = paste("~/plink-1.07-x86_64/plink", "--bfile", first_file, "--merge-list", written_file, "--make-bed", "--out", out_file, sep = " "))
misssnp_file <- file.path(current_dir, "data_integrated", "23andme_TR-merge.missnp")
for (file_name in complete_valid_binary_files) {
  
  in_file <- file.path(current_dir, "data_binary", file_name)
  out_file <- file.path(current_dir, "data_binary_temp", file_name)
  system(command = paste("~/plink-1.07-x86_64/plink", "--bfile", in_file, "--exclude", misssnp_file, "--allow-no-sex", "--make-bed", "--out", out_file, sep = " "))
  
}
bed_files <- file.path(current_dir, "data_binary_temp", paste0(complete_valid_binary_files, ".bed"))
bim_files <- file.path(current_dir, "data_binary_temp", paste0(complete_valid_binary_files, ".bim"))
fam_files <- file.path(current_dir, "data_binary_temp", paste0(complete_valid_binary_files, ".fam"))
df_filenames <- data.frame(bed = bed_files, bim = bim_files, fam = fam_files)
first_file <- file.path(current_dir, "data_binary_temp", complete_valid_binary_files[1])
df_filenames <- df_filenames[c(-1), ]
write.table(x = df_filenames, file = file.path(current_dir, "data_integrated", "merge_files.txt"), sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)
system(command = paste("~/plink-1.07-x86_64/plink", "--bfile", first_file, "--merge-list", written_file, "--make-bed", "--allow-no-sex", "--out", out_file, sep = " "))
