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




# Docking result clustering
def get_RMSD(ref, target):
    distances = []
    maxSubStr = rdFMCS.FindMCS([ref, target])
    a = ref.GetSubstructMatch(Chem.MolFromSmarts(maxSubStr.smartsString))
    b = target.GetSubstructMatch(Chem.MolFromSmarts(maxSubStr.smartsString))   
    for atomA,atomB in list(zip(a, b)):
        pos_A = ref.GetConformer().GetAtomPosition(atomA)
        pos_B = target.GetConformer().GetAtomPosition(atomB)
        coord_A = np.array((pos_A.x, pos_A.y, pos_A.z))
        coord_B = np.array((pos_B.x, pos_B.y, pos_B.z))
        dist_numpy = np.linalg.norm(coord_A - coord_B)        
        distances.append(dist_numpy)
    rmsd = np.round(math.sqrt(1/len(distances)*sum([i*i for i in distances])),3)
    return rmsd

def get_score(inputL, inputE, result="full"):
    ID = os.path.basename(inputL)[:-11]
    df = PandasTools.LoadSDF(inputL)
    df["ID"] = [ID + "_" + df.loc[i,"POSE"] for i in df.index]
    df["RMSD_EL"] = [get_RMSD(inputE, df.loc[i, "ROMol"]) for i in df.index]
    df["SCORE"] = pd.to_numeric(df["SCORE"])
    df["RMSD_EL"] = pd.to_numeric(df["RMSD_EL"])
    df["RMSD_LB_BP"] = pd.to_numeric(df["RMSD_LB_BP"])
    df["RMSD_UB_BP"] = pd.to_numeric(df["RMSD_UB_BP"])
    first = ["ID", "SCORE", "RMSD_EL"]
    full = ["ID", "POSE", "SCORE", "RMSD_EL", "RMSD_LB_BP", "RMSD_UB_BP"]
    if result == "first": 
        return df[first].iloc[:1]
    if result == "full": 
        return df[full]
