# Install necessary packages if not already installed

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite")
}
if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table
  ")
}

# Load required libraries
library(jsonlite)
library(data.table)

time <- function(...) { 
  time_measurement <- system.time(eval(...)) 
  time_measurement[["user.self"]] 
} 
# -------------------------
# CAS to CID
# -------------------------


t <- time(pc_mapping <- fread("/storage/skall/CAS_to_SMILES/CAS_data_wrangler/PubChem/CID-CAS-filtered", header=FALSE,sep=" ", col.names = c("CID","CAS") ))

pc_mapping$Synonyms <- trimws(pc_mapping$CAS)

# Remove duplicates
pc_mapping <- unique(pc_mapping)

# Create a lookup table: a list where each CID points to its CAS number(s)
lookup_table <- split(pc_mapping$CID, pc_mapping$CAS)

json_file_path <- "/storage/skall/CAS_to_SMILES/CAS_data_wrangler/PubChem/CAS-CID-lookup_table.json"
t <- time(write_json(lookup_table, path = json_file_path, pretty = TRUE))
print("time", t)
