#!/usr/bin/env python
# coding: utf-8

# In[5]:
from ppc_froi_project.Reaction_Functions_python import *
from ppc_froi_project.Lists_python import *
from ppc_froi_project.Dictionnaries_python import *


from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from IPython.display import display
import sys

from rdkit.Chem import PandasTools
#import pandas as pd
from pathlib import Path
import os
from rdkit.Chem.Draw import IPythonConsole
IPythonConsole.ipython_useSVG=True

#get_ipython().run_line_magic('run', './Reaction_Functions.ipynb')

#import pytest 

############# function needed for the test  #############

def canonicalize_smi(smi):
    return Chem.MolToSmiles(Chem.MolFromSmiles(smi))

############# Reaction function tests #############
def test_carbonyles_deprotonation_function():
    start_smiles = "c1ccccc1CC(=O)C1CCCCC1"
    result = Chem.MolToSmiles(carbonyles_deprotonation_function(start_smiles),canonical=True)
    assert result == canonicalize_smi("[O-]C(=Cc1ccccc1)C1CCCCC1")

def test_Imine_Formation2_function():
    start_smiles = "C1ccC(=O)CC1"
    result = Chem.MolToSmiles(Imine_Formation2_function(start_smiles),canonical=True)
    assert result == canonicalize_smi("C1ccC(=N)CC1")

def test_addition_cyanide_function ():
    start_smiles = "c1ccccc1CC(=O)C1CCCCC1"
    result = Chem.MolToSmiles(addition_cyanide_function(start_smiles),canonical=True)
    assert result == canonicalize_smi("N#CC(O)(Cc1ccccc1)C1CCCCC1")

def test_heterocycle_synthesis():
    start_smiles = "C(C)(C)C(=O)C(C)CC(=O)c1ccccc1"
    start_smiles2 = "CN"
    result = Chem.MolToSmiles(heterocycle_synthesis_function(start_smiles, start_smiles2),canonical=True)
    assert result == canonicalize_smi("Cc1cc(-c2ccccc2)n(C)c1C(C)C")

def test_hydrolysis_cyanohydrine_function():
    start_smiles = "N#CC(O)(Cc1ccccc1)C1CCCCC1"
    result = Chem.MolToSmiles(hydrolysis_cyanohydrine_function(start_smiles),canonical=True)
    assert result == canonicalize_smi("O=C(O)C(O)(Cc1ccccc1)C1CCCCC1")

def test_OxineFormation_function():
    start_smiles = "c1C(=O)cCc1"
    start_smiles2 = "NO"
    result = Chem.MolToSmiles(OxineFormation_function(start_smiles, start_smiles2),canonical=True)
    assert result == canonicalize_smi("c1C(=NO)cCc1")

def test_wolfkischner_function():
    start_smiles = "C1ccCC(=O)CCC1"
    result = Chem.MolToSmiles(WolfKishner_function(start_smiles), canonical=True)
    assert result == canonicalize_smi("C1ccCCCCC1")

def strecker_reaction_function(input_smiles):
    start_smiles = "c1ccccc1CC(=O)CCC(=O)CCC"
    result = Chem.MolToSmiles(strecker_reaction_function(start_smiles),canonical=True)
    assert result == canonicalize_smi("c1ccccc1CC(O)(N)CCC(=O)CCC")
    
def pinacol_function (input_smiles):
    start_smiles = "OC1(C2(O)CCCC2)CCCC1"
    result = Chem.MolToSmiles(pinacol_function (start_smiles),canonical=True)
    assert result == canonicalize_smi("C2CCCC12C(=O)CCCC1")
    
def Hydration_function(input_smiles):
    start_smiles = "C=[O+]CCC=O"
    result = Chem.MolToSmiles(Hydration_function(start_smiles),canonical=True)
    assert result == canonicalize_smi("C=[O+]CCC(O)(O)")
    
def bisulfite_add_compounds1_function(input_smiles):
    start_smiles = "c1ccccc1CC(=O)CCC(=O)CCC"
    result = Chem.MolToSmiles(bisulfite_add_compounds1_function(start_smiles),canonical=True)
    assert result == canonicalize_smi("c1ccccc1CC(O)(S(=O)(=O)[O-].[Na+])CCC(=O)CCC")
    
def reduction_cyanohydrine_function(input_smiles):
    start_smiles = "c1ccccc1CC(=O)CCC(=O)CCC"
    result = Chem.MolToSmiles(reduction_cyanohydrine_function(start_smiles),canonical=True)
    assert result == canonicalize_smi("c1ccccc1CC(O)(CN)CCC(=O)CCC")
    
############# Test for function needed in some reaction #############

def test_neighbor_index():
    start_smiles = "NC(=O)CC"
    atom = "O"
    start_mol = Chem.MolFromSmiles(start_smiles)
    neighbor = neighbor_index(start_mol, atom)
    assert neighbor == 3

