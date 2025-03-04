# script.R
library(optparse)
library(dplyr)
library(tidyr)
library(data.table)
library(pbapply) 
pbo = pboptions(type="txt")

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="tmp_out.txt", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("--cas_cid_file"), type="character", default="CID-CAS-filtered",
              help="pubchem cid cas mapping [default= %default]", metavar="character"),
  make_option(c("--cid_smiles_file"), type="character", default="CID-SMILES",
              help="pubchem cid smiles mapping [default= %default]", metavar="character"),
  make_option(c("--debug"), type="logical", default=FALSE,
              help="debug with prints and 10 rows [default= %default]", metavar="logical")
) 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# SCRIPT STARTS HERE ############################################
 
# -------------------------
# CAS to CID
# -------------------------
start_time <- Sys.time()

# Step 1: Read the CSV Dataframe
print("Loading data...")

if (grepl("\\.zip$", opt$file)) {
  # Extract CSV filename from ZIP (assumes it has the same name but with .csv extension)
  csv_file <- sub("\\.csv.zip$", ".csv", basename(opt$file))
  
  df <- read.csv(
    unz(opt$file, csv_file), 
    header = TRUE, 
    sep = ",", 
    quote = "\"", 
    fill = TRUE, 
    stringsAsFactors = FALSE, 
    na.strings = c("", "NA")
  )
} else {
  df <- read.csv(
    opt$file, 
    header = TRUE, 
    sep = ",", 
    quote = "\"", 
    fill = TRUE, 
    stringsAsFactors = FALSE, 
    na.strings = c("", "NA")
  )
}

# Step 2: Split rows with multiple CAS numbers into duplicate rows
if (opt$debug){
  df <- head(df, 10)
  df <- df %>%
  mutate(CAS = strsplit(as.character(CAS), ";")) %>% # Split CAS numbers by ';'
  unnest(CAS) %>% # Unnest the list into separate rows
  mutate(CAS = trimws(CAS)) 
  print(head(df))
} else {
  df <- df %>%
  mutate(CAS = strsplit(as.character(CAS), ";")) %>% # Split CAS numbers by ';'
  unnest(CAS) %>% # Unnest the list into separate rows
  mutate(CAS = trimws(CAS)) # Trim any whitespace around CAS numbers
}

# Load data with CID and CAS
print("Loading CID --> CAS mapping...")
pc_mapping <- fread(opt$cas_cid_file, header=FALSE,sep=" ", col.names = c("CID","CAS"))
if (opt$debug){
  print(head(pc_mapping))
}
setkey(pc_mapping, CAS)

# Function to map CAS to CID
map_cas_to_cid <- function(cas) {
  tryCatch({
    cid <- pc_mapping[.(cas), nomatch = NA]$CID[1]
    return(cid)
  }, error = function(e) {
    return(NA)
  })
}

# Step 2: Get unique CAS numbers
unique_cas <- unique(df$CAS)

# Add a CID column by mapping CAS to CID
print(paste("Fetching", length(unique_cas), "CIDs..."))
cid_mapping <- pbsapply(unique_cas, map_cas_to_cid, USE.NAMES = TRUE)
# Filter out NA values from cid_mapping
valid_cid_mapping <- cid_mapping[!is.na(cid_mapping)]

print(paste("Failed to find CID for", length(unique_cas)-length(valid_cid_mapping), "CAS numbers..."))

# Convert the result to a data frame for easier joining
cid_df <- data.frame(CAS = unique_cas, CID = cid_mapping, row.names=NULL, stringsAsFactors = FALSE)
df <- df %>%
  left_join(cid_df, by = "CAS")

print("Loading CID --> SMILES mapping...")
if (opt$debug){
  pc_mapping <- fread(opt$cid_smiles_file, header=FALSE, sep="\t", col.names = c("CID","SMILES"), nrows=100000)
  print(head(pc_mapping))
} else {
  pc_mapping <- fread(opt$cid_smiles_file, header=FALSE, sep="\t", col.names = c("CID","SMILES"))
}
setkey(pc_mapping, CID)

# Function to map CID to SMILES
map_cid_to_smiles <- function(cid) {
  tryCatch({
    smiles <- pc_mapping[.(cid), nomatch = NA]$SMILES[1]
    return(smiles)
  }, error=function(e){
    print(e)
    return(NA)
    })
}

# Step 5: Add a SMILES column by mapping CID to SMILES
unique_cid <- unique(df$CID)

# Add a SMILES column by mapping CID to SMILES
print(paste("Fetching", length(unique_cid), "SMILES..."))
smiles_mapping <- pbsapply(unique_cid, map_cid_to_smiles, USE.NAMES = TRUE)
# Filter out NA values from smiles_mapping
valid_smiles_mapping <- smiles_mapping[!is.na(smiles_mapping)]
print(paste("Failed to find SMILES for", length(unique_cid)-length(valid_smiles_mapping), "CID numbers..."))

# Convert the result to a data frame for easier joining
smiles_df <- data.frame(CID = unique_cid, SMILES_PC = smiles_mapping, row.names=NULL, stringsAsFactors = FALSE)
df <- df %>%
  left_join(smiles_df, by = "CID")

# Step 6: Output the final dataframe
if (opt$debug){
  print(head(df[, c("CAS", "CID", "SMILES_PC")]))
}

# Optionally save the final dataframe to a new CSV file
if (grepl("\\.zip$", opt$out)) {
  # Define CSV filename (same as opt$out but with .csv extension)
  csv_file <- sub("\\.csv.zip$", ".csv", opt$out)
  
  # Write the CSV file
  write.csv(df, csv_file, row.names = FALSE)
  
  # Zip the CSV file
  zip::zipr(opt$out, files = csv_file)
  
  # Optionally, remove the original CSV after zipping
  file.remove(csv_file)
} else {
  # Write as usual if not zipping
  write.csv(df, opt$out, row.names = FALSE)
}

execution_time <- Sys.time() - start_time
print(paste("Execution time:", execution_time))