import os
import sys


import ast
import math
import plip
import timeit
import shutil
import py3Dmol
import contextlib
import xlsxwriter
import urllib.request

import numpy as np
import pandas as pd

from google.colab import drive, files
from tqdm.notebook import tqdm
from openbabel import pybel
from Bio.PDB import PDBIO, PDBParser
from rdkit import Chem
from rdkit.Chem import rdFMCS, AllChem, PandasTools
from plip.exchange.report import BindingSiteReport
from plip.structure.preparation import PDBComplex




# Binding interactions profiler
def interaction_profiler(inputPL):
    protlig = PDBComplex()
    INTER = {}
    INTER_TYPE = []
    ID = os.path.basename(inputPL[:-7])
    CSV_Ifile = os.path.join(INT_FLD, f"{ID[:-2]}/{ID}_inter.csv")
    COL = ["RESNR", "RESTYPE", "RESCHAIN", "RESNR_LIG", "RESTYPE_LIG", 
           "RESCHAIN_LIG", "DIST", "LIGCARBONIDX", "PROTCARBONIDX", "LIGCOO", 
           "PROTCOO"]
    INTER_DF = pd.DataFrame(columns=COL)

    protlig.load_pdb(inputPL)
    LIG = protlig.ligands[0]
    protlig.characterize_complex(LIG)

    BD_keys = BOND_DICT.keys()
    BSI_key = f":".join([str(x) for x in LIG.members[0]])
    BSI_val = protlig.interaction_sets[BSI_key]
    BSR = BindingSiteReport(BSI_val) 
    INTER[BSI_key] = {k: [getattr(BSR, f"{k}_features")] + 
                      getattr(BSR, f"{k}_info")
                      for k in tuple(BD_keys)}
    for BD_key in list(BD_keys):
        DF = pd.DataFrame.from_records(
            INTER[BSI_key][BD_key][1:],
            columns=INTER[BSI_key][BD_key][0])
        if not DF.empty:
            INTER_TYPE.extend([BD_key.upper()]*len(DF))
            INTER_DF = INTER_DF.append(DF)
    INTER_DF = INTER_DF.assign(BOND=INTER_TYPE)
    INTER_DF.reset_index(drop=True, inplace=True)
    INTER_DF.to_csv(CSV_Ifile, index=False)

def view_interaction(inputCSV, result="summary"):
    DIST_CALC = []
    COLOR = []
    MIDCOO = []
    interaction = pd.read_csv(inputCSV, converters={
        "LIGCOO": lambda x: ast.literal_eval(str(x)), 
        "PROTCOO": lambda x: ast.literal_eval(str(x)), 
        "BOND": lambda x: x.lower()})
    for LC, PC, BT in zip(interaction["LIGCOO"], interaction["PROTCOO"], 
                          interaction["BOND"]):
        COLOR.append(BOND_DICT[BT][1])
        p1 = np.array([LC[0], LC[1], LC[2]])
        p2 = np.array([PC[0], PC[1], PC[2]])
        squared_dist = np.sum((p1-p2)**2, axis=0)
        dist = np.round(np.sqrt(squared_dist) ,2)
        DIST_CALC.append(dist)
        mid_x = np.round((LC[0] + PC[0]) / 2, 2)
        mid_y = np.round((LC[1] + PC[1]) / 2, 2)
        mid_z = np.round((LC[2] + PC[2]) / 2, 2)
        p_mid = (mid_x, mid_y, mid_z) 
        MIDCOO.append(p_mid)
    interaction["BOND"] = interaction["BOND"].str.upper()
    interaction["COLOR"] = COLOR
    interaction["MIDCOO"] = MIDCOO
    interaction["DIST_CALC"] = DIST_CALC
    interaction["RESNAME"] = interaction["RESTYPE"] + interaction["RESNR"].astype(str)

    summary = ["RESNR", "RESTYPE", "DIST_CALC", "BOND", "COLOR"]
    py3dmol = ["RESNR", "DIST_CALC", "LIGCOO", "PROTCOO", "MIDCOO", "BOND"]
    if result == "summary":
        df = interaction[summary]
    if result == "py3Dmol":
        df = interaction[py3dmol]
    return df
