import os
import sys
import numpy as np
import pandas as pd
import math
from Bio.PDB import PDBIO, PDBParser
from rdkit import Chem
from rdkit.Chem import rdFMCS, AllChem, PandasTools


class Hide:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, "w")  
    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

# Biological assembly detector and extractor
def chain_profilier(inputFile):
    parser = PDBParser()
    structure = parser.get_structure("X", inputFile)
    chain_obj = [s for s in structure.get_chains()] 
    chain_IDs = [c.get_id() for c in chain_obj]
    chain_len = len(chain_obj)
    return [chain_len, chain_IDs, chain_obj]

def seperate_chain(inputFile):
    io = PDBIO()
    name = os.path.basename(inputFile)[:-4]
    dir = os.path.dirname(inputFile)
    chain_info = chain_profilier(inputFile)
    chain_len, chain_IDs, chain_obj = chain_info[0], chain_info[1], chain_info[2]
    print(f"+ Chains detected: {chain_len} ({', '.join(chain_IDs)})")
    for ID,chain in zip(chain_IDs, chain_obj):
        chain_pdb = name + "_" + ID + ".pdb"
        chain_pdb_file = os.path.join(dir, chain_pdb)
        io.set_structure(chain)
        io.save(chain_pdb_file)
        print(f"+ {chain_pdb} ==> {os.path.basename(dir)} folder")
