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

def test_strecker_reaction_function():
    start_smiles = "c1ccccc1CC(=O)C1CCCCC1"
    result = Chem.MolToSmiles(strecker_reaction_function(start_smiles),canonical=True)
    assert result == canonicalize_smi("c1ccccc1CC(O)(N)C1CCCCC1")

def test_pinacol_function ():
    start_smiles = "C1CCC1(O)C(O)(c1ccccc1)(c1ccccc1)"
    result = Chem.MolToSmiles(pinacol_function (start_smiles),canonical=True)
    assert result == canonicalize_smi("C1CCC(=O)C1(c1ccccc1)(c1ccccc1)")

def test_pinacol_function1 ():
    start_smiles = "C1CCCC1(O)C(O)(c1ccccc1)(c1ccccc1)"
    result = Chem.MolToSmiles(pinacol_function (start_smiles),canonical=True)
    assert result == canonicalize_smi("C1CCCC(=O)C1(c1ccccc1)(c1ccccc1)")

def test_pinacol_function2 ():
    start_smiles = "C1CCCCC1(O)C(O)(c1ccccc1)(c1ccccc1)"
    result = Chem.MolToSmiles(pinacol_function (start_smiles),canonical=True)
    assert result == canonicalize_smi("C1CCCCC(=O)C1(c1ccccc1)(c1ccccc1)")

def test_pinacol_function3 ():
    start_smiles = "C1CCCCCC1(O)C(O)(c1ccccc1)(c1ccccc1)"
    result = Chem.MolToSmiles(pinacol_function (start_smiles),canonical=True)
    assert result == canonicalize_smi("C1CCCCCC(=O)C1(c1ccccc1)(c1ccccc1)")

def test_pinacol_function4 ():
    start_smiles = "C1CCCCCCC1(O)C(O)(c1ccccc1)(c1ccccc1)"
    result = Chem.MolToSmiles(pinacol_function (start_smiles),canonical=True)
    assert result == canonicalize_smi("C1CCCCCCC(=O)C1(c1ccccc1)(c1ccccc1)")

def test_pinacol_function5 ():
    start_smiles = "C1CC1(O)C(O)(c1ccccc1)(c1ccccc1)"
    result = Chem.MolToSmiles(pinacol_function (start_smiles),canonical=True)
    assert result == canonicalize_smi("C1CC(=O)C1(c1ccccc1)(c1ccccc1)")

def test_Hydration_function():
    start_smiles = "C1CCC(=O)CC1"
    result = Chem.MolToSmiles(Hydration_function(start_smiles),canonical=True)
    assert result == canonicalize_smi("C1CCC(O)(O)CC1")
    
def test_bisulfite_add_compounds1_function():
    start_smiles = "c1ccccc1CC(=O)C1CCCCC1"
    result = Chem.MolToSmiles(bisulfite_add_compounds1_function(start_smiles),canonical=True)
    assert result == canonicalize_smi("c1ccccc1CC(O)(S(=O)(=O)[O-].[Na+])C1CCCCC1")
    
def test_reduction_cyanohydrine_function():
    start_smiles = "c1ccccc1CC(=O)C1CCCCC1"
    result = Chem.MolToSmiles(reduction_cyanohydrine_function(start_smiles),canonical=True)
    assert result == canonicalize_smi("c1ccccc1CC(O)(CN)C1CCCCC1")
    
############# Test for function needed in some reaction #############

def test_neighbor_index():
    start_smiles = "NC(=O)CC"
    atom = "O"
    start_mol = Chem.MolFromSmiles(start_smiles)
    neighbor = neighbor_index(start_mol, atom)
    assert neighbor == 3

