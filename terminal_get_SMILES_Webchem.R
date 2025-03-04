# script.R
library(optparse)
library(dplyr)
library(data.table)
library(pbapply)  
library(webchem)
pbo = pboptions(type="txt")

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="tmp_out.txt", 
              help="output file name [default= %default]", metavar="character"),
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
# CAS to SMILES
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

if (opt$debug){
  df <- head(df, 10)
  print(df[, c("CAS","CID","SMILES_PC")])
}

# Function to map CAS to CID
map_cas_to_smiles <- function(cas) {
  tryCatch({
    smiles <- cir_query(identifier=cas, representation="smiles", match="first")[2]
    return(smiles)
  }, error = function(e) {
    return(NA)
  })
}

# Step 1: Filter out rows where CAS is NA or "-"
filtered_df <- df %>%
  filter(!is.na(CAS) & CAS != "-" & is.na(SMILES_PC))

# Step 2: Get unique CAS numbers
unique_cas <- unique(filtered_df$CAS)

# Step 3: Apply the function to unique CAS numbers
# Use sapply to apply the function to each unique CAS number
print(paste("Fetching", length(unique_cas), "SMILES..."))
smiles_mapping <- pbsapply(unique_cas, map_cas_to_smiles, USE.NAMES = TRUE)

# Convert the result to a data frame for easier joining
smiles_df <- data.frame(CAS = unique_cas, SMILES_CIR = unlist(smiles_mapping), row.names=NULL, stringsAsFactors = FALSE)

# Step 4: Merge the SMILES back into the original dataframe
df <- df %>%
  left_join(smiles_df, by = "CAS")

# Step 6: Output the final dataframe
# Print the first few rows of the final dataframe
if (opt$debug){
  print(head(df[, c("CAS", "CID", "SMILES_PC", "SMILES_CIR")]))
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