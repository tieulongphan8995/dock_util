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
