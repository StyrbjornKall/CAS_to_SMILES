# Install necessary packages if not already installed

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite")
}

# Load required libraries
library(jsonlite)

# -------------------------
# CID to SMILES
# -------------------------

# Read the Pubchem file
pc_mapping <- read.table(
  "/storage/skall/CAS_to_SMILES/CAS_data_wrangler/PubChem/CID-SMILES", 
  header=FALSE, 
  sep="\t", 
  col.names = c("CID","SMILES")
)

pc_mapping$SMILES <- trimws(pc_mapping$SMILES)

# Remove duplicates
pc_mapping <- unique(pc_mapping)

# Create a lookup table: a list where each CID points to its CAS number(s)
lookup_table <- split(pc_mapping$CID, pc_mapping$SMILES)

json_file_path <- "/storage/skall/CAS_to_SMILES/CAS_data_wrangler/PubChem/CID-SMILES-lookup_table.json"
write_json(lookup_table, path = json_file_path, pretty = TRUE)