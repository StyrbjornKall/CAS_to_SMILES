#!/bin/bash

# Ensure necessary commands are available
for cmd in awk grep; do
  if ! command -v $cmd &> /dev/null; then
    echo "Error: $cmd is not installed. Please install it."
    exit 1
  fi
done

# Variables
input_file="/storage/skall/CAS_to_SMILES/CAS_data_wrangler/PubChem/CID-Synonym-filtered"
output_file="/storage/skall/CAS_to_SMILES/CAS_data_wrangler/PubChem/CID-CAS-filtered"
temp_file=$(mktemp)
cas_pattern="^[0-9]{2,7}-[0-9]{2}-[0-9]$"

# Count the total number of lines in the input file
total_lines=$(wc -l < "$input_file")
current_line=0

# Process the file and provide a progress indication
awk -v pattern="$cas_pattern" -v total_lines="$total_lines" '
BEGIN { FS = "\t" }
{
  if (match($2, pattern)) {
    print $1, $2
  }
}
' "$input_file" | tee "$temp_file"

# Sort and remove duplicates (if necessary)
sort -u "$temp_file" > "$output_file"

# Clean up
rm "$temp_file"

# Final progress output
echo "Extraction complete. Results written to $output_file"