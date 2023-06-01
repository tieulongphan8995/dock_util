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

# Gridbox methods
def get_atom_coordinate(inputE):
    mol = Chem.MolFromPDBFile(inputE)
    molAtom = mol.GetAtoms()
    molConf = mol.GetConformer()
    xcoor = [molConf.GetAtomPosition(c).x for c in range(len(molAtom))]
    ycoor = [molConf.GetAtomPosition(c).y for c in range(len(molAtom))]
    zcoor = [molConf.GetAtomPosition(c).z for c in range(len(molAtom))]
    return [xcoor, ycoor, zcoor]

def labogrid(coorList, scale=2):
    X, Y, Z = coorList[0], coorList[1], coorList[2]
    range = [[min(X), max(X)], 
             [min(Y), max(Y)], 
             [min(Z), max(Z)]]
    center = [np.round(np.mean(range[0]), 3), 
              np.round(np.mean(range[1]), 3), 
              np.round(np.mean(range[2]), 3)]
    size = [np.round(abs(range[0][0]-range[0][1])*scale, 3),
            np.round(abs(range[1][0]-range[1][1])*scale, 3),
            np.round(abs(range[2][0]-range[2][1])*scale, 3)]
    return [center, size]

def eboxsize(coorList, scale=1):
    gyBoxRatio = 0.23
    X, Y, Z = coorList[0], coorList[1], coorList[2]
    center = [np.round(np.mean(X), 3),
              np.round(np.mean(Y), 3),
              np.round(np.mean(Z), 3),]
    dist = [(x-center[0])**2 + 
            (y-center[1])**2 + 
            (z-center[2])**2 
            for x,y,z in zip(X,Y,Z)]
    size = [np.round((math.sqrt(np.sum(dist)/len(X))/gyBoxRatio)*scale, 3)]*3
    return [center, size]

def eboxsize_mod(coorList, scale=1):
    gyBoxRatio = 0.23
    X, Y, Z = coorList[0], coorList[1], coorList[2]
    range = [[min(X), max(X)], 
             [min(Y), max(Y)], 
             [min(Z), max(Z)]]
    center = [np.round(np.mean(range[0]), 3), 
              np.round(np.mean(range[1]), 3), 
              np.round(np.mean(range[2]), 3)]
    dist = [(x-center[0])**2 + 
            (y-center[1])**2 + 
            (z-center[2])**2 
            for x,y,z in zip(X,Y,Z)]
    size = [np.round((math.sqrt(np.sum(dist)/len(X))/gyBoxRatio)*scale, 3)]*3
    return [center, size]

def autodock_grid(coorList, scale=1):
    min_len = 22.5
    X, Y, Z = coorList[0], coorList[1], coorList[2]
    range = [[min(X), max(X)], 
             [min(Y), max(Y)], 
             [min(Z), max(Z)]]
    center = [np.round(np.mean(range[0]), 3), 
              np.round(np.mean(range[1]), 3), 
              np.round(np.mean(range[2]), 3)]
    size = [min_len]*3
    return [center, size]

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

# Py3Dmol model viewer
def VIEW_PROT(inputP, model="Cartoon", color="white", focusRes="",
              focusResColor="white", addStick=False, addLine=False,
              showChain=False, showResLabel=False, showVDW=False, 
              outline=False):
    # Variables
    mview = py3Dmol.view(1000,1500)
    if model in "none":
        model = cType = color = ""
    if model in "cartoon":
        cType = "color"
    if model in ("line", "stick"):
        cType, color = "colorscheme", color + "Carbon"

    # Protein Model
    count = 1
    prot_model = count
    print(f"+ Showing {os.path.basename(inputP)}")
    mol = open(inputP, "r").read()
    mview.addModel(mol, "pdb")
    mview.setStyle(
        {"model": prot_model},
        {model: {cType: color}})
    
    # Show chains
    chainLen = chain_profilier(inputP)[0]
    if showChain and chainLen > 1:
        for n,c in zip(range(chainLen), COLORS):
            mview.setStyle(
                {"and": [{"model": prot_model}, {"chain": chr(65+n)}]},
                {model: {cType: c if model == "cartoon" else c + "Carbon"}})
            mview.addLabel(
                f"Chain {chr(65+n)}",
                {"fontColor": c, "backgroundOpacity": 0.7, "alignment": "topLeft"},
                {"and": [{"model": prot_model}, {"chain": chr(65+n)}]})
  
    # Focus RES
    if focusRes != "":
        res = focusRes.split(",")
        mview.addStyle(
            {"and": [{"model": prot_model}, {"resi": res}]},
            {"stick": {"colorscheme": focusResColor + "Carbon", "radius": 0.15}})
        # Label focused RES
        if showResLabel:
            mview.addResLabels(
                {"and": [{"model": prot_model}, {"resi": res}]},
                {"alignment": "bottomLeft",
                 "showBackground": False,
                 "inFront": True,
                 "fixed": False,
                 "fontSize": 14,
                 "fontColor": "0x000000",
                 "screenOffset": {"x": 15, "y": 15}})
  
    # Show VDW surface
    if showVDW:
        mview.addSurface(
            py3Dmol.VDW, 
            {"color": "white", "opacity": 0.4},
            {"model": prot_model})

    # Show outline
    if outline:
        mview.setViewStyle(
            {"style": "outline", 
            "color": "black", 
            "width": 0.1})
  
    print(f"")
    mview.setBackgroundColor("0xFFFFFF")
    mview.center({"model": prot_model})
    mview.show()

def VIEW_ELIG(inputE, showAtomLabel=False, outline=False):
    # Variable
    mview = py3Dmol.view(1000, 1500)

    # Experimental ligand model
    count = 1
    elig_model = count
    print(f"+ Showing {os.path.basename(inputE)}")
    mol = open(inputE, "r").read()
    mview.addModel(mol, "pdb")
    mview.setStyle(
        {"model": elig_model},
        {"stick": {"colorscheme": "lightGrayCarbon"}})
  
    # Label all atoms
    if showAtomLabel:
        mview.addPropertyLabels(
            "atom",
            {"model": elig_model},
            {"backgroundOpacity": 0.7, "inFront": True})
    
    # Show outline
    if outline:
        mview.setViewStyle(
            {"style": "outline", 
            "color": "black", 
            "width": 0.1})

    print(f"")
    mview.setBackgroundColor("0xFFFFFF")
    mview.zoomTo({"model": elig_model})
    mview.show()

def VIEW_LIG(inputL, showHs=False, showAtomLabel=False, outline=False):
    # Variable
    mview = py3Dmol.view(1000, 1500)

    # Ligand model
    count = 1
    lig_model = count
    print(f"+ Showing {os.path.basename(inputL)}")
    smi = open(inputL, "r").read()
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)

    # Remove Hs
    if not showHs:
        mol = Chem.RemoveHs(mol)
    mblock = Chem.MolToMolBlock(mol)
    mview.addModel(mblock, "mol")
    mview.setStyle(
        {"model": lig_model},
        {"stick": {"colorscheme": "lightGrayCarbon"}})
  
    # Label all atoms
    if showAtomLabel:
        mview.addPropertyLabels(
            "atom",
            {"model": lig_model},
            {"backgroundOpacity": 0.7, "inFront": True})
  
    # Show outline
    if outline:
        mview.setViewStyle(
            {"style": "outline", 
            "color": "black", 
            "width": 0.1})

    print(f"")
    mview.setBackgroundColor("0xFFFFFF")
    mview.center({"model": lig_model})
    mview.show()

def VIEW_GRID(inputP, inputE, focusRes, center, size=[10, 10, 10]):
    # Variable
    mview = py3Dmol.view(1000, 1500)

    # Grid box
    bxi, byi, bzi = center[0], center[1], center[2]
    bxf, byf, bzf = size[0], size[1], size[2]
    print(f"+ Center: X {center[0]}  Y {center[1]}  Z {center[2]}")
    print(f"+ Size: W {size[0]}  H {size[1]}  D {size[2]}")
    mview.addBox(
        {"center":{"x":bxi, "y":byi, "z":bzi}, 
        "dimensions": {"w":bxf, "h":byf, "d":bzf}, 
        "color": "skyBlue", "opacity": 0.7})

    # Protein model
    count = 1
    prot_model = count
    print(f"+ Showing {os.path.basename(inputP)}")
    molA = open(inputP, "r").read()
    mview.addModel(molA, "pdb")
    mview.setStyle(
        {"model": prot_model},
        {"cartoon": {"color": "white"}})
    
    # Experimental ligand model
    count += 1
    elig_model = count
    print(f"+ Showing {os.path.basename(inputE)}")
    molB = open(inputE, "r").read()
    mview.addModel(molB, "pdb")
    mview.setStyle(
        {"model": elig_model},
        {"stick": {"colorscheme": "greenCarbon"}})
  
    # Focus RES
    if focusRes != "":
        res = focusRes.split(",")
        mview.addStyle(
            {"and": [{"model": prot_model}, {"resi": res}]}, 
            {"stick": {"colorscheme": "orangeCarbon", "radius": 0.15}})
        mview.addResLabels(
            {"and": [{"model": prot_model}, {"resi": res}]},
            {"alignment": "bottomLeft",
            "showBackground": False,
            "inFront": True,
            "fixed": False,
            "fontSize": 14,
            "fontColor": "0x000000",
            "screenOffset": {"x": 15, "y": 15}})

    print(f"")
    mview.setBackgroundColor("0xFFFFFF")
    mview.center({"model": elig_model})
    mview.show()

def VIEW_PILE(inputP, inputPP, inputL, inputE, inputCSV, model="Cartoon", 
              color="white", focusRes="", focusResColor="white", 
              showInter=False, showSubunit=False, showResLabel=False, 
              showVDW=False, showPartProt=False, showExpLig=False, 
              showLig=False, showAllPose=False, showBestPose=False,
              outline=False):
    # Variables
    resUQ = []
    mview = py3Dmol.view(1000,1500)
    if model in "none":
        model = cType = color = ""
    if model in "cartoon":
        cType = "color"
    if model in ("line", "stick"):
        cType, color = "colorscheme", color + "Carbon"

    # Protein model
    count = 1 
    prot_model = count
    print(f"+ Showing {os.path.basename(inputP)}")
    molA = open(inputP, "r").read()
    mview.addModel(molA, "pdb")
    mview.setStyle(
        {"model":prot_model},
        {model: {cType: color}})

  # Show chains
    chainLen = chain_profilier(inputP)[0]
    if showSubunit and chainLen > 1:
        for n, c in zip(range(chainLen), COLORS):
            mview.setStyle(
                {"and": [{"model": prot_model}, {"chain": chr(65+n)}]},
                {model: {cType: c if model == "cartoon" else c + "Carboon"}})
            mview.addLabel(
                f"Subunit {chr(65+n)}",
                {"fontColor": c, "backgroundOpacity": 0.7, "alignment": "topLeft"},
                {"and": [{"model": prot_model}, {"chain": chr(65+n)}]})
            
    # Focus RES
    if focusRes != "":
        resUQ.extend([int(r) for r in focusRes.split(",")])
        mview.addStyle(
            {"and": [{"model": prot_model}, {"resi": resUQ}]},
            {"stick": {"colorscheme": focusResColor + "Carbon", "radius": 0.15}})
        # Label focused RES
        if showResLabel:
            mview.addResLabels(
                {"and": [{"model": prot_model}, {"resi": resUQ}]},
                {"alignment": "bottomLeft",
                 "showBackground": False,
                 "inFront": True,
                 "fixed": False,
                 "fontSize": 14,
                 "fontColor": "0x000000",
                 "screenOffset": {"x": 15, "y": 15}})
    
    # Interactions
    if showInter:
        interaction = view_interaction(inputCSV, result="py3Dmol")
        RESNR = interaction["RESNR"]
        DIST_CALC = interaction["DIST_CALC"]
        LIGCOO = interaction["LIGCOO"]
        PROTCOO = interaction["PROTCOO"]
        MIDCOO = interaction["MIDCOO"]
        BOND = interaction["BOND"]
        for RN,DC,LC,PC,MC,BT in zip(RESNR, DIST_CALC, LIGCOO, PROTCOO, MIDCOO, BOND):
            BT = BT.lower()
            if RN not in resUQ:
                resUQ.append(RN)
                mview.addStyle(
                    {"and": [{"model": prot_model}, {"resi": RN}]},
                    {"stick": {"colorscheme": "whiteCarbon", "radius": 0.15}})
                mview.addResLabels(
                    {"and": [{"model": prot_model}, {"resi": RN}]},
                    {"alignment": "bottomLeft", "showBackground": False,
                     "inFront": True, "fixed": False, "fontSize": 14, 
                     "fontColor": "0x000000", "screenOffset": {"x": 15, "y": 15}})
                mview.addCylinder(
                    {"start": {"x": LC[0], "y": LC[1], "z": LC[2]},
                     "end": {"x": PC[0], "y": PC[1], "z": PC[2]}, 
                     "radius": 0.05, "fromCap": 1, "toCap": 1, 
                     "color": BOND_DICT[BT][0], "dashed": True})
                mview.addLabel(
                    str(DC) + " Å",
                    {"position": {"x": MC[0], "y": MC[1], "z": MC[2]},
                     "alignment": "bottomLeft", "inFront": False, "fixed": False,
                     "backgroundColor": BOND_DICT[BT][0], "fontSize": 10,
                     "screenOffset": {"x": 5, "y": 5}})

    # Show VDW surface
    if showVDW:
        mview.addSurface(
            py3Dmol.VDW, 
            {"color": "white", "opacity": 0.4},
            {"model": prot_model})

    # Partner protein model
    if showPartProt and os.path.exists(inputPP):
        count += 1
        pprot_model = count
        print(f"+ Showing {os.path.basename(inputPP)} (yellow)")
        molB = open(inputPP, "r").read()
        mview.addModel(molB, "pdb")
        mview.setStyle(
            {"model": pprot_model},
            {"cartoon": {"color": "yellow"}})
        mview.addStyle(
            {"model": pprot_model},
            {"stick": {"colorscheme": "yellowCarbon"}})
  
    # Experimental ligand model
    if showExpLig:
        count += 1
        elig_model = count
        print(f"+ Showing {os.path.basename(inputE)} (gray)")
        molC = open(inputE, "r").read()
        mview.addModel(molC, "pdb")
        mview.setStyle(
            {"model": elig_model},
            {"stick": {"color": "gray"}})
    
    # Ligand model
    if showLig:
        count += 1
        lig_model = count
        print(f"+ Showing {os.path.basename(inputL)} (orange)")
        molD = open(inputL, "r").read()
        mview.addModel(molD, "pdb")
        mview.setStyle(
            {"model": lig_model},
            {"stick": {"colorscheme": "orangeCarbon"}})

    # Show poses of selected ligand 
    if showAllPose:
        pose = sorted([ os.path.join(os.path.dirname(inputL), f) for f in os.listdir(os.path.dirname(inputL)) if f.endswith(".pdb") ])
        pose.remove(inputL[len(DCK_FLD)+1:]) if inputL[len(DCK_FLD)+1:] in pose else None
        for f in pose:
            count += 1
            pose_model = count
            fs = "".join(f)
            mol5 = open(os.path.join(DCK_FLD, fs), "r").read()
            mview.addModel(mol5, "pdb")
            mview.setStyle(
                {"model": pose_model},
                {"stick": {"color": "blue", "opacity": 0.5, "radius": 0.2}})

    if outline:
        mview.setViewStyle(
            {"style": "outline", "color": "black", "width": 0.1})
    
    print(f"")
    mview.setBackgroundColor("0xFFFFFF")
    mview.center({"model": lig_model})
    mview.enableFog(True)
    mview.show()
