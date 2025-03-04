# Map CAS to SMILES
This repo contains a bunch of convenience scripts for mapping a list of CAS numbers to SMILES. It does some preprocessing and makes a lookup using a local copy of PubChem, and API calls through Webchem (Cactus) and PubMed.

Based on ~100 000 CAS numbers tested, PubChem usually finds 95% of them. About 1000 of the remaining ones could be found in Webchem and another 200 through PubMed.