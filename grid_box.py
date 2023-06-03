import os
import sys
import math
import timeit
import shutil
import contextlib
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFMCS, AllChem, PandasTools
import py3Dmol



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

  while True:
      try:
        from pymol import cmd
        def getbox(selection='sele', extending = 6.0):

          ([minX, minY, minZ],[maxX, maxY, maxZ]) = cmd.get_extent(selection)

          minX = minX - float(extending)
          minY = minY - float(extending)
          minZ = minZ - float(extending)
          maxX = maxX + float(extending)
          maxY = maxY + float(extending)
          maxZ = maxZ + float(extending)

          SizeX = maxX - minX
          SizeY = maxY - minY
          SizeZ = maxZ - minZ
          CenterX =  (maxX + minX)/2
          CenterY =  (maxY + minY)/2
          CenterZ =  (maxZ + minZ)/2

          cmd.delete('all')

          return {'center_x':CenterX,'center_y': CenterY, 'center_z': CenterZ},{'size_x':SizeX,'size_y': SizeY,'size_z': SizeZ}
        break
      except:
        print("No module Pymol")
        break
