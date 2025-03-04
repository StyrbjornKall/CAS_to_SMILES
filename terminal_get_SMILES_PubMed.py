import requests
import pandas as pd
import argparse
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
import time
import random

# Create a session for better performance
session = requests.Session()

def cas_to_smiles_pubmed(cas, debug=False):
    
    """Fetch SMILES from PubChem API with error handling and retries."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{cas}/property/IsomericSMILES/JSON"
    
    for attempt in range(3):  # Retry up to 3 times
        try:
            response = session.get(url, timeout=5)  # Set timeout
            response.raise_for_status()  # Raise error for bad response
            data = response.json()
            smiles = data['PropertyTable']['Properties'][0]['IsomericSMILES']
            return cas, smiles
        except requests.exceptions.HTTPError as e:
            if response.status_code == 404:
                if debug:
                    print(f"[PubChem] CAS {cas} not found.")
                return cas, None
        except requests.exceptions.HTTPError as e:
            if response.status_code == 404:
                if debug:
                    print(f"[PubChem] CAS {cas} not found.")
                return cas, None
        except requests.exceptions.RequestException as e:
            if debug:
                print(f"[PubChem] Network error for CAS {cas}: {e}")
        time.sleep(random.uniform(1, 3))  # Random sleep to avoid hitting API too fast
    return cas, None  # Return None after retries
    
def cas_to_smiles_cactus(cas, debug=False):
    """Fetch SMILES from Cactus API with error handling."""
    url = f"https://cactus.nci.nih.gov/chemical/structure/{cas}/smiles"

    try:
        response = session.get(url, timeout=5)
        response.raise_for_status()
        smiles = response.text.strip()
        return cas, smiles if smiles else None
    except requests.exceptions.RequestException as e:
        if debug:
            print(f"[Cactus] Error fetching CAS {cas}: {e}")
        return cas, None

def fetch_smiles_parallel(cas_list, func, debug=False, max_workers=4):
    """Fetch SMILES in parallel using ThreadPoolExecutor."""
    results = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_cas = {executor.submit(func, cas, debug): cas for cas in cas_list}
        for future in tqdm(as_completed(future_to_cas), total=len(cas_list), desc="Fetching SMILES"):
            cas, smiles = future.result()
            results.append((cas, smiles))
    return results

def test():
    cas_numbers = [
        "50-00-0",  # Formaldehyde
        "67-64-1",  # Acetone
        "71-43-2",  # Benzene
        "64-17-5",  # Ethanol
        "7732-18-5"  # Water
    ]
    
    molecules = [cas_to_smiles_pubmed(cas) for cas in cas_numbers ] 
    
    df = pd.DataFrame()
    df['CAS'] = cas_numbers
    df['SMILES'] = molecules
    
    output_file = 'cas_to_smiles_results.csv'
    df.to_csv(output_file, index=False)
    
    print("\nResults:")
    print(df)

def main():
    parser = argparse.ArgumentParser(description='Process and pickle a file.')
    parser.add_argument('--file', required=True, help='Input file path')
    parser.add_argument('--out', required=True, help='Output file path')
    parser.add_argument('--debug', required=False, default=False, help='debug less data')
    parser.add_argument('--cir', required=False, default=False, help='do cir query as well')

    args = parser.parse_args()

    if (args.debug == 'TRUE') or (args.debug == 'True') or (args.debug == 'true'):
        debug = True
    else:
        debug = False

    if (args.cir == 'TRUE') or (args.cir == 'True') or (args.cir == 'true'):
        cir = True
    else:
        cir = False
    
    # Read the input file
    if args.file.endswith('.zip'):
        df = pd.read_csv(args.file, compression='zip', low_memory=False)
    else:
        df = pd.read_csv(args.file, low_memory=False)
    
    if debug:
        df = df.head(100)

    unique_cas = df[~df.CAS.isna() & df.SMILES_PC.isna()]['CAS'].unique()
    
    print(f'Number of unique CAS numbers with no SMILES: {len(unique_cas)}/{df.CAS.nunique()}')

    print(f"Fetching {len(unique_cas)} SMILES from PubMed...")
    smiles_mapping_pm = fetch_smiles_parallel(unique_cas, cas_to_smiles_pubmed, debug=debug)
    # Convert to DataFrame
    smiles_df_pm = pd.DataFrame(smiles_mapping_pm, columns=['CAS', 'SMILES_PM'])
    df = df.merge(smiles_df_pm, on='CAS', how='left')

    print(f'Failed to find SMILES for {len(unique_cas)-smiles_df_pm.SMILES_PM.dropna().nunique()} CAS numbers.')
    unique_cas = df[~df.CAS.isna() & df.SMILES_PC.isna() & df.SMILES_PM.isna()]['CAS'].unique()
    print(f'Number of unique CAS numbers with no SMILES: {len(unique_cas)}/{df.CAS.nunique()}')

    if cir:
        unique_cas_missing = df[df['SMILES_PM'].isna()]['CAS'].dropna().unique()
        
        print(f"Fetching {len(unique_cas)} SMILES from Cactus...")
        smiles_mapping_cactus = fetch_smiles_parallel(unique_cas_missing, cas_to_smiles_cactus, debug=args.debug)

        smiles_df_cactus = pd.DataFrame(smiles_mapping_cactus, columns=['CAS', 'SMILES_CACTUS'])
        df = df.merge(smiles_df_cactus, on='CAS', how='left')
        
        print(f'Failed to find SMILES for {len(unique_cas)-smiles_df_cactus.SMILES_CACTUS.dropna().nunique()} CAS numbers.')

    if debug:
        print(df.head())

    # Save the output file
    if args.out.endswith('.zip'):
        df.to_csv(args.out, compression='zip', index=False)
    else:
        df.to_csv(args.out, index=False)

if __name__ == "__main__":
    main()
