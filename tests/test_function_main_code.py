#!/usr/bin/env python
# coding: utf-8

# In[5]:

from ppc_froi_project.Report import *
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
#import pytest

############# function needed for the test  #############

def canonicalize_smi(smi):
    return Chem.MolToSmiles(Chem.MolFromSmiles(smi))

############# Function of the code #############

def test_validate_smiles():
    smiles = "COCOCO"
    return_smiles = validate_smiles(smiles)
    correct_smiles = canonicalize_smi(smiles)
    assert Chem.MolToSmiles(return_smiles) == correct_smiles

########################################### EXEMPLES A RAJOUTER 3333333333333333333333333333333333333333333333333333333333333333333

def test_subgroup_in_start():
    mol = Chem.MolFromSmiles("C(=O)C(=O)C")
    subgroup_list = subgroup_in_start(mol, Start)
    correct_list = ['CC(=O)C', 'CC(=O)']
    assert subgroup_list == correct_list

def test_subgroup_in_finish():
    mol = Chem.MolFromSmiles("C(O)C(O)C")
    subgroup_list = subgroup_in_finish(mol, Finish)
    correct_list = ['CC(O)C', 'CC(O)C', 'CCO', 'CCC', 'CO']
    assert subgroup_list == correct_list

def test_check_for_cond():
    condition = "Acid"
    check_condition = check_for_cond(condition, Condition)
    assert check_condition == True

def test_check_same_dict_in_list():
    start = "CC(=O)C"
    condition = "Acid and Nucleophile"
    start_mol = Chem.MolFromSmiles(start)
    name_test = check_same_dict_in_list(start_mol, condition, Reaction)
    correct_name = 'nucleophilic_attack2'
    assert name_test == correct_name

def test_check_same_dict_in_list1():
    start = Chem.MolFromSmiles("CC(=O)C")
    start2 = Chem.MolFromSmiles("CN")
    name_test = check_same_dict_in_list1(start, start2, Reaction)
    correct_name = 'Imine_Formation1'
    assert name_test == correct_name


def check_same_dict_in_list3():
    finish = Chem.MolFromSmiles("CC(O)C")
    condition = "Base and Nucleophile"
    name_test = check_same_dict_in_list3(finish, condition, Reaction)
    correct_name = 'nucleophilic_attack1'
    assert name_test == correct_name

#def test_check_same_dict_in_list4():
 #   condition =
  #  finish =
   # finish2 =
    #name_test = check_same_dict_in_list4(condition, finish, finish2, Reaction)
    #assert name_test == correct_name

def test_check_same_dict_in_list5():
    start = Chem.MolFromSmiles("CC(=O)CCC(=O)C")
    start2 = Chem.MolFromSmiles("CN")
    condition = "C2OOH"
    name_test = check_same_dict_in_list5(start,start2,condition,Reaction)
    correct_name = 'heterocycle_synthesis'
    assert name_test == correct_name


# In[ ]:




