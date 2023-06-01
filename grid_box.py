import os
import sys
import math
import timeit
import shutil
import contextlib
import xlsxwriter
import urllib.request
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFMCS, AllChem, PandasTools




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
