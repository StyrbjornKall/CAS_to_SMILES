import pandas as pd
import argparse
import json

def main():
    parser = argparse.ArgumentParser(description='Process and pickle a file.')
    parser.add_argument('--file', required=True, help='Input file path')
    parser.add_argument('--out', required=True, help='Output file path')
    parser.add_argument('--debug', required=False, default=False, help='debug less data')
    parser.add_argument('--drop_missing_smiles', required=False, default=False, help='drops rows with missing SMILES')
    parser.add_argument('--save_cas_smiles_dict', required=False, default=True, help='saves a lookup json file for CAS to SMILES')

    args = parser.parse_args()

    if (args.debug == 'TRUE') or (args.debug == 'True') or (args.debug == 'true'):
        debug = True
    else:
        debug = False

    if (args.save_cas_smiles_dict == 'TRUE') or (args.save_cas_smiles_dict == 'True') or (args.save_cas_smiles_dict == 'true'):
        save_cas_smiles_dict = True
    else:
        save_cas_smiles_dict = False
    
    if (args.drop_missing_smiles == 'TRUE') or (args.drop_missing_smiles == 'True') or (args.drop_missing_smiles == 'true'):
        drop_missing_smiles = True
    else:
        drop_missing_smiles = False
    
    # Read the input file
    if args.file.endswith('.zip'):
        df = pd.read_csv(args.file, compression='zip', low_memory=False)
    else:
        df = pd.read_csv(args.file, low_memory=False)
    
    if debug:
        df = df.head(100)

    print('Joining SMILES...')
    df['SMILES'] = df['SMILES_PC']
    if 'SMILES_CIR' in df.columns:
        df['SMILES'] = df['SMILES'].fillna(df['SMILES_CIR'])
        df.drop(columns=['SMILES_PC', 'SMILES_CIR'], inplace=True)
    if 'SMILES_PM' in df.columns:
        df['SMILES'] = df['SMILES'].fillna(df['SMILES_PM'])
        df.drop(columns=['SMILES_PM'], inplace=True)
    if 'SMILES_CACTUS' in df.columns:
        df['SMILES'] = df['SMILES'].fillna(df['SMILES_CACTUS'])
        df.drop(columns=['SMILES_CACTUS'], inplace=True)

    if drop_missing_smiles:
        df.dropna(subset=['SMILES'], inplace=True)

    # Save the output file
    if args.out.endswith('.zip'):
        df.to_csv(args.out, compression='zip', index=False)
    else:
        df.to_csv(args.out, index=False)

    if save_cas_smiles_dict:
        cas_smiles_dict = df[['CAS', 'SMILES']].set_index('CAS').to_dict()['SMILES']
        if args.out.endswith('.zip'):
            with open(args.out.replace('.csv.zip', '_cas_smiles_dict.json'), 'w') as f:
                json.dump(cas_smiles_dict, f)
        else:
            with open(args.out.replace('.csv', '_cas_smiles_dict.json'), 'w') as f:
                json.dump(cas_smiles_dict, f)

if __name__ == "__main__":
    main()
