#!/bin/bash
# Specify the interpreter
source /local/Anaconda2-2.4.1-Linux-x86_64/etc/profile.d/conda.sh

conda activate r_env_4.2
# Define paths to files
input_file="/storage/skall/tox_data/preprocessed/Preprocessed_tox_data_SK_20250212.csv.zip"
pubchem_cas_cid_mapping="/storage/skall/CAS_to_SMILES/CAS_data_wrangler/PubChem/CID-CAS-filtered"
pubchem_cid_smiles_mapping="/storage/skall/CAS_to_SMILES/CAS_data_wrangler/PubChem/CID-SMILES"
output_file="/storage/skall/tox_data/preprocessed/Preprocessed_tox_data_SK_20250212_with_SMILES.csv.zip"

debug=FALSE

# Get CIDs and SMILES from Pubchem
Rscript terminal_get_CID_SMILES_PubChem.R --file=$input_file --out="tmp_1.csv" --cas_cid_file=$pubchem_cas_cid_mapping --cid_smiles_file=$pubchem_cid_smiles_mapping --debug=$debug

# Get missing from webchem
Rscript terminal_get_SMILES_Webchem.R --file="tmp_1.csv" --out="tmp_2.csv" --debug=$debug

conda deactivate
conda activate trident_311

# Get all SMILES with pubmed and cactus manually
python terminal_get_SMILES.py --file="tmp_2.csv" --out="tmp_3.csv" --debug=$debug --cir=FALSE
python terminal_concatenate_smiles.py --file="tmp_3.csv" --out=$output_file --debug=$debug --drop_missing_smiles=FALSE

# Save file as pkl in python
python /home/skall/Preprocess_toxicity_data/pickle_files.py --file=$output_file

conda deactivate
