#!/usr/bin/env python
# coding: utf-8

# <div style="font-size:24px;"><b> Reaction Functions </b></div>

# <div style="font-size:11pt;"> 
# This section is dedicated to showcasing our functions that are designed to compute and return the final product. These functions have been developed using a variety of Python libraries, demonstrating diverse programming approaches to achieve the desired outcomes.</div>

# <span style="font-size:14pt; text-decoration:underline;"> Necessary Imports </span>

# In[3]:


from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import PandasTools
from rdkit.Chem.Draw import IPythonConsole
#import pandas as pd
from pathlib import Path
import os


# <span style="font-size:14pt; text-decoration:underline;"> Complementary Functions </span>

# In[17]:

## Find the index of the second attached carbon on finish_pattern ##
def neighbor_index(finish,atom):
    oxygene_idx = int(finish.GetSubstructMatches(Chem.MolFromSmiles(atom))[0][0])                         
    neighbor_carbon_idx = finish.GetAtomWithIdx(oxygene_idx).GetNeighbors() 
    carbon_neighbor_idx = None
    
    # Get the neihgbor's neighbor of the oxygene atom #
    for neighbor in neighbor_carbon_idx:
        neighbors_neighbor = neighbor.GetNeighbors()
    # Get idx of the neighbor's neighbor of the O atom #
        for neighbor_neighbor in neighbors_neighbor:
            if neighbor_neighbor.GetIdx() != oxygene_idx:
                carbon_neighbor_idx = neighbor_neighbor.GetIdx()
    return carbon_neighbor_idx

# <span style="font-size:14pt; text-decoration:underline;"> Deprotonation of Carbonyl </span>

# In[4]:


def carbonyles_deprotonation_function(input_smiles):

############################# Take the molec apart #########################################

    ## set the variable ##
    start_mol = Chem.MolFromSmiles(input_smiles)
    start_pattern = Chem.MolFromSmiles("C=O")
    finish_pattern = Chem.MolFromSmiles("C[O-]")

    ## slicing the molecules to find R1 and R2 ##
    R12 = Chem.DeleteSubstructs(start_mol, start_pattern)
    R12_smiles = Chem.MolToSmiles(R12)
    # slicing for R1 and R2 #
    R_list = R12_smiles.split(".")
    R1_smiles = R_list[0]
    R1 = Chem.MolFromSmiles(R1_smiles)

    ## find R1, R2 and start idx ##
    match_start_mol_pattern = start_mol.GetSubstructMatch(start_pattern)
    R1_idx = 0
    R2_idx = len(R1.GetAtoms())
    finish_idx = len(R12.GetAtoms()) 

############################# Molecule ... ASSEMBLE! #########################################

    ### String Assemble ##
    big_smiles = R12_smiles + "." + Chem.MolToSmiles(finish_pattern)
    big_smiles_mol = Chem.MolFromSmiles(big_smiles)
    big_smiles_mol_editable = Chem.EditableMol(big_smiles_mol)
    #map_big_smiles = mol_with_atom_and_molecule_index(big_smiles_mol)
    
    ## Making bond ##
    big_smiles_mol_editable.AddBond(R1_idx, finish_idx, order=Chem.rdchem.BondType.SINGLE)
    big_smiles_mol_editable.AddBond(R2_idx, finish_idx, order=Chem.rdchem.BondType.DOUBLE)
    finish_mol1 = big_smiles_mol_editable.GetMol()

    ## Sanitize ##
    Chem.SanitizeMol(finish_mol1)

#########################################################################################

    return finish_mol1


# <span style="font-size:14pt; text-decoration:underline;"> Imine Formation under Acidic Conditions  </span>

# In[5]:


def Imine_Formation2_function(start_input_smiles):

##Reaction Information   
    Imine_Formation2 = { "Name": "Addition of amines on acidic conditions",
        "Start": "C(=O)C", 
        "Condition": "4<pH<6",
        "Finish": "C(=N)C"}  
    
## Input Transformation from smiles to molecule   
    input_mol = Chem.MolFromSmiles(start_input_smiles)
    start_mol = Chem.MolFromSmiles(Imine_Formation2["Start"])

## Pattern Identification  
    if input_mol.HasSubstructMatch(start_mol):
        for match in input_mol.GetSubstructMatches(start_mol):
            edit_mol = Chem.EditableMol(input_mol)
            _, carbonyl_idx, _ = match  # Match indices: start C, C=O, end C
            
## Molecule Assemble 
            edit_mol.ReplaceAtom(carbonyl_idx, Chem.Atom(7))
            modified_mol = edit_mol.GetMol()
            
## Sanitize
            Chem.SanitizeMol(modified_mol)

##Output 
            final_mol = Chem.MolToSmiles(modified_mol, isomericSmiles=True)
            molecule = Chem.MolFromSmiles(final_mol)
            
            return molecule 


# <span style="font-size:14pt; text-decoration:underline;"> Imine Formation under Basic Conditions </span>

# In[6]:


def Imine_Formation1_function(start_input_smiles):

##Reaction Information
    Imine_Formation1 = { "Name": "Addition of amines on basic conditions" ,
                        "Start": "CC(=O)C" , 
                        "Start2":"NC",
                        "Condition" : "pH>6", 
       
                        "Finish":"CC(=NC)C" }
    
## Input Transformation from smiles to molecule  
    input_mol = Chem.MolFromSmiles(start_input_smiles)
    start_mol = Chem.MolFromSmiles(Imine_Formation2["Start"])

## Pattern Identification
    if input_mol.HasSubstructMatch(start_mol):
        for match in input_mol.GetSubstructMatches(start_mol):
            edit_mol = Chem.EditableMol(input_mol)
            _, carbonyl_idx, _ = match  # Match indices: start C, C=O, end C

## Molecule Assemble           
            edit_mol.ReplaceAtom(carbonyl_idx, Chem.Atom(7))
            modified_mol = edit_mol.GetMol()
            
##Sanitize          
            Chem.SanitizeMol(modified_mol)
            
## Output           
            final_mol = Chem.MolToSmiles(modified_mol, isomericSmiles=True)
            molecule = Chem.MolFromSmiles(final_mol)
            
            return molecule 


# <span style="font-size:14pt; text-decoration:underline;"> Addition of Cyanide </span>

# In[7]:


def addition_cyanide_function (input_smiles):
    ############################# Take the molec apart #########################################

    ## set the variable ##
    start_mol = Chem.MolFromSmiles(input_smiles)
    start_pattern = Chem.MolFromSmiles("C=O")
    finish_pattern = Chem.MolFromSmiles("C(O)(C#N)")

    ## slicing the molecules to find R1 and R2 ##
    R12 = Chem.DeleteSubstructs(start_mol, start_pattern)
    R12_smiles = Chem.MolToSmiles(R12)
    # slicing for R1 and R2 #
    R_list = R12_smiles.split(".")
    R1_smiles = R_list[0]
    R1 = Chem.MolFromSmiles(R1_smiles)

    ## find R1, R2 and start idx ##
    match_start_mol_pattern = start_mol.GetSubstructMatch(start_pattern)
    R1_idx = 0
    R2_idx = len(R1.GetAtoms())
    finish_idx = len(R12.GetAtoms()) + 2

############################# Molecule ... ASSEMBLE! #########################################

    ### String Assemble ##
    big_smiles = R12_smiles + "." + Chem.MolToSmiles(finish_pattern)
    big_smiles_mol = Chem.MolFromSmiles(big_smiles)
    big_smiles_mol_editable = Chem.EditableMol(big_smiles_mol)
    #map_big_smiles = mol_with_atom_and_molecule_index(big_smiles_mol)
    
    ## Making bond ##
    big_smiles_mol_editable.AddBond(R1_idx, finish_idx, order=Chem.rdchem.BondType.SINGLE)
    big_smiles_mol_editable.AddBond(R2_idx, finish_idx, order=Chem.rdchem.BondType.SINGLE)
    finish_mol1 = big_smiles_mol_editable.GetMol()
    
    ## Sanitize ##
    Chem.SanitizeMol(finish_mol1)
    finish_mol_smiles = Chem.MolToSmiles(finish_mol1)

#########################################################################################

    return finish_mol1


# <span style="font-size:14pt; text-decoration:underline;"> Hydrolysis of Cyanohydrine </span>

# In[8]:


def hydrolysis_cyanohydrine_function(input_smiles):
    ############################# Take the molec apart #########################################

    ## set the variable ##
    start_mol = Chem.MolFromSmiles(input_smiles)
    start_pattern = Chem.MolFromSmiles("C(O)(C#N)")
    finish_pattern = Chem.MolFromSmiles("C(O)(C(=O)O)")
    
    ## slicing the molecules to find R1 and R2 ##
    R12 = Chem.DeleteSubstructs(start_mol, start_pattern)
    R12_smiles = Chem.MolToSmiles(R12)

    # slicing for R1 and R2 #
    R_list = R12_smiles.split(".")
    R1_smiles = R_list[0]
    R1 = Chem.MolFromSmiles(R1_smiles)
    R2_smiles = R_list[1]
    R2 = Chem.MolFromSmiles(R2_smiles)

    ## find R1, R2 and start idx ##
    match_start_mol_pattern = start_mol.GetSubstructMatch(start_pattern)
    R1_idx = 0
    R2_idx = len(R1.GetAtoms())
    finish_idx = len(R12.GetAtoms()) + len(finish_pattern.GetAtoms()) - 2

############################# Molecule ... ASSEMBLE! #########################################

    ### String Assemble ##
    big_smiles = R12_smiles + "." + Chem.MolToSmiles(finish_pattern)
    big_smiles_mol = Chem.MolFromSmiles(big_smiles)
    big_smiles_mol_editable = Chem.EditableMol(big_smiles_mol)
    #map_big_smiles = mol_with_atom_and_molecule_index(big_smiles_mol)
    
    ## Making bond ##
    big_smiles_mol_editable.AddBond(R1_idx, finish_idx, order=Chem.rdchem.BondType.SINGLE)
    big_smiles_mol_editable.AddBond(R2_idx, finish_idx, order=Chem.rdchem.BondType.SINGLE)
    finish_mol1 = big_smiles_mol_editable.GetMol()
    finish_mol_smiles = Chem.MolToSmiles(finish_mol1)
    ## Sanitize ##
    #Chem.SanitizeMol(finish_mol1)
    
    finish_mol_smiles = Chem.MolToSmiles(finish_mol1)
    
#########################################################################################
    
    return finish_mol1


# <span style="font-size:14pt; text-decoration:underline;"> Formation of Oxine </span>

# In[9]:


def OxineFormation_function(start_input_smiles, start2_input_smiles):

## Reaction Information
    oxine_formation = {
        "Name": "Addition of nucleophilic nitrogen to form an oxine",
        "Start": "C(=O)",  
        "Start2": "NO",
        "Finish": "C(=NO)" 
    }
## Input Transformation from smiles to molecule 
    molecule1 = Chem.MolFromSmiles(start_input_smiles)
    molecule2 = Chem.MolFromSmiles(start2_input_smiles)

## Pattern Identification 
    pattern1 = Chem.MolFromSmarts(oxine_formation["Start"])  
    nitroso_pattern = Chem.MolFromSmarts(oxine_formation["Start2"]) 
    
    if molecule1.HasSubstructMatch(pattern1) and molecule2.HasSubstructMatch(nitroso_pattern):
        
## Molecule Assemble
        modified_smiles = start_input_smiles.replace("C(=O)", "C(=NO)")
        
## Output
        final_mol = Chem.MolFromSmiles(modified_smiles)
        
        return final_mol


# <span style="font-size:14pt; text-decoration:underline;"> Wolf-Kishner-Huang Carbonyl Reduction  </span>

# In[10]:


def WolfKishner_function(start_input_smiles):

##Reaction Information
    wolf_kishner_huang_reduction = { "Finish" : "CCC",
                                    "Condition" : "NH2NH2, KOH, 180°C",
                                    "Start" : "CC(=O)C", 
                                    "Name" : "Wolf-Kishner-Huang Reduction" }
    molecule = Chem.MolFromSmiles(start_input_smiles)

##Pattern Identification
    pattern = Chem.MolFromSmarts(wolf_kishner_huang_reduction["Start"])
    if molecule.HasSubstructMatch(pattern):
        
##Molecule Assemble 
        finish_smiles = start_input_smiles.replace("C(=O)", "C")
        
##Output
        final_mol = Chem.MolFromSmiles(finish_smiles)
    
        return final_mol


# <span style="font-size:14pt; text-decoration:underline;"> Strecker Reaction  </span>

# In[19]:


def strecker_reaction_function(input_smiles):
     ############################# Take the molec apart #########################################

    ## set the variable ##
    start_mol = Chem.MolFromSmiles(input_smiles)
    start_pattern = Chem.MolFromSmiles("C=O")
    finish_pattern = Chem.MolFromSmiles("C(O)(N)")
    
    ## slicing the molecules to find R1 and R2 ##
    R12 = Chem.DeleteSubstructs(start_mol, start_pattern)
    R12_smiles = Chem.MolToSmiles(R12)

    # slicing for R1 and R2 #
    R_list = R12_smiles.split(".")
    R1_smiles = R_list[0]
    R1 = Chem.MolFromSmiles(R1_smiles)
    R2_smiles = R_list[1]
    R2 = Chem.MolFromSmiles(R2_smiles)

    ## find R1, R2 and start idx ##
    match_start_mol_pattern = start_mol.GetSubstructMatch(start_pattern)
    R1_idx = 0
    R2_idx = len(R1.GetAtoms())
    finish_idx = len(R12.GetAtoms()) + len(finish_pattern.GetAtoms()) - 2

    ############################# Molecules ... ASSEMBLE! #########################################

    ### String Assemble ##
    big_smiles = R12_smiles + "." + Chem.MolToSmiles(finish_pattern)
    big_smiles_mol = Chem.MolFromSmiles(big_smiles)
    big_smiles_mol_editable = Chem.EditableMol(big_smiles_mol)
    
    ## Making bond ##
    big_smiles_mol_editable.AddBond(R1_idx, finish_idx, order=Chem.rdchem.BondType.SINGLE)
    big_smiles_mol_editable.AddBond(R2_idx, finish_idx, order=Chem.rdchem.BondType.SINGLE)
    finish_mol1 = big_smiles_mol_editable.GetMol()
    
    ## Sanitize ##
    Chem.SanitizeMol(finish_mol1)
    
    finish_mol_smiles = Chem.MolToSmiles(finish_mol1) 
    
    return finish_mol1


# <span style="font-size:14pt; text-decoration:underline;"> Pinacol Rearrangement  </span>

# In[12]:


def pinacol_function (input_smiles):
############################# Take the molec apart #########################################

    ## set the variable ##
    start_mol = Chem.MolFromSmiles(input_smiles)

    ## select how much carbon there is in the cycle ##
    if len(start_mol.GetSubstructMatch(Chem.MolFromSmiles("C1CC1(O)C(O)"))) == 6:
        start_pattern = Chem.MolFromSmiles("C1CC1(O)C(O)")
        finish_pattern = Chem.MolFromSmiles("C1CCC1(=O)")
        
    ## select how much carbon there is in the cycle ##
    if len(start_mol.GetSubstructMatch(Chem.MolFromSmiles("C1CCC1(O)C(O)"))) == 7:
        start_pattern = Chem.MolFromSmiles("C1CCC1(O)C(O)")
        finish_pattern = Chem.MolFromSmiles("C1CCCC1(=O)")
        
    elif len(start_mol.GetSubstructMatch(Chem.MolFromSmiles("C1CCCC1(O)C(O)"))) == 8:
        start_pattern = Chem.MolFromSmiles("C1CCCC1(O)C(O)")
        finish_pattern = Chem.MolFromSmiles("C1CCCCC1(=O)")

    elif len(start_mol.GetSubstructMatch(Chem.MolFromSmiles("C1CCCCC1(O)C(O)"))) == 9:
        start_pattern = Chem.MolFromSmiles("C1CCCCC1(O)C(O)")
        finish_pattern = Chem.MolFromSmiles("C1CCCCCC1(=O)")
        
    elif len(start_mol.GetSubstructMatch(Chem.MolFromSmiles("C1CCCCCC1(O)C(O)"))) == 10:
        start_pattern = Chem.MolFromSmiles("C1CCCCCC1(O)C(O)")
        finish_pattern = Chem.MolFromSmiles("C1CCCCCCC1(=O)")
        
    elif len(start_mol.GetSubstructMatch(Chem.MolFromSmiles("C1CCCCCCC1(O)C(O)"))) == 11:
        start_pattern = Chem.MolFromSmiles("C1CCCCCCC1(O)C(O)")
        finish_pattern = Chem.MolFromSmiles("C1CCCCCCCC1(=O)")
        
    ## Extract R's from the start molec ##
    R12 = Chem.DeleteSubstructs(start_mol, start_pattern)
    R12_smiles = Chem.MolToSmiles(R12)
    # chekc if there is a R2 #
    if "." in R12_smiles:
        R12_list = R12_smiles.split(".")
        R1_smiles = R12_list[0]
        R2_smiles = R12_list[1]
    else:
        R1_smiles = R12_smiles
    
    ### Find the index of the R's attached atoms ##
    big_smiles = R12_smiles + "." + Chem.MolToSmiles(finish_pattern)
    R1_idx = 0
    
    # check if there is a R2 #
    if "." in R12_smiles:
        R2_idx = len(Chem.MolFromSmiles(R2_smiles).GetAtoms()) + R1_idx 
    
    ## Find the index of the attached carbon on finish_pattern ##            
    carbon_idx = neighbor_index(Chem.MolFromSmiles(big_smiles),"O")
    
############################# Molecules ... ASSEMBLE! #########################################    

    ## Modify big_smiles to become editable ##
    big_smiles_mol = Chem.MolFromSmiles(big_smiles)
    big_smiles_mol_editable = Chem.EditableMol(big_smiles_mol)
    
    ## create a bond between R1 and finish_pattern ##
    big_smiles_mol_editable.AddBond(R1_idx, carbon_idx, order=Chem.rdchem.BondType.SINGLE)
    
    # check if there is a R2 #
    if "." in R12_smiles:
        
        big_mol = big_smiles_mol_editable.GetMol()
        Chem.SanitizeMol(big_mol)
        carbon_idx2 = neighbor_index(big_mol,"O")
        big_smiles_mol_editable2 = Chem.EditableMol(big_mol)
        big_smiles_mol_editable2.AddBond(R2_idx, carbon_idx2, order=Chem.rdchem.BondType.SINGLE)
        finish_mol1 = big_smiles_mol_editable2.GetMol()
        
    else:    
        finish_mol1 = big_smiles_mol_editable.GetMol()

    ## Last Treatment ##
    Chem.SanitizeMol(finish_mol1)
    finish_mol_smiles = Chem.MolToSmiles(finish_mol1)
    
#########################################################################################

    return finish_mol1


# <span style="font-size:14pt; text-decoration:underline;"> Hydration of Carbonyl </span>

# In[13]:


def Hydration_function(input_smiles):

##Reaction Information
    hydration = {
        "Name": "Hydration of Carbonyls",
        "Start": "C=O",   
        "Condition": "H2O and acid",
        "Finish": "[O-][H+]" 
    }
##Pattern Identification
    molecule = Chem.MolFromSmiles(input_smiles)
    pattern = Chem.MolFromSmarts(hydration["Start"])  
    
    if molecule.HasSubstructMatch(pattern):
##Molecule Assemble 
        finish_mol = Chem.MolFromSmiles("C(O)(O)")
        modified_mols = AllChem.ReplaceSubstructs(molecule, pattern, finish_mol, replaceAll=True)
        final_mol = Chem.MolToSmiles(modified_mols[0], isomericSmiles=True)
## Output 
        molecule = Chem.MolFromSmiles(final_mol)
        
        return molecule 


# <span style="font-size:14pt; text-decoration:underline;"> Heterocycle Synthesis : Pyrrole Formation </span>

# In[14]:


def heterocycle_synthesis_function(start, start2):
    
    ## set the variables ##
    start_pattern = Chem.MolFromSmiles("C(=O)CCC(=O)")
    start_mol = Chem.MolFromSmiles(start)
    finish_pattern = Chem.MolFromSmiles("c1cc[N]c1")
    start2_mol = Chem.MolFromSmiles(start2)
    start2_pattern = Chem.MolFromSmiles("N")
    
############################# Finding Index of atoms which reacts #########################################
   
    ## Get the big smiles and mol by adding start and start2 together ##
    big_smiles = start + "." + start2
    big_mol = Chem.MolFromSmiles(big_smiles)
    
    
    ## Get the index of the oxygen's neighbor carbon ##
    # Get the list of the oxygen, index in start molecule #
    oxygen_tuple = big_mol.GetSubstructMatches(Chem.MolFromSmiles("O"))
    oxygen_list = []
    for i in range(len(oxygen_tuple)):
        oxygen_list.append(oxygen_tuple[i][0])
    number_oxygen = len(oxygen_list)
    
    # Find the list of the two oxygen's index corresponding to the start pattern #
    if number_oxygen > 2:
        pattern_tuple = big_mol.GetSubstructMatches(start_pattern)
        pattern_list = []
        for i in range(len(pattern_tuple[0])):
            pattern_list.append(pattern_tuple[0][i])
        pattern_set = set(pattern_list)
        oxygen_set = set(oxygen_list)
        list_idx = list(pattern_set.intersection(oxygen_set))
        oxygen_list = list_idx
    
    # Finally get the index of carbons #
    list_carbon_idx = []
    for i  in range(len(oxygen_list)):
        #neighbor_carbon_idx = big_mol.GetAtomWithIdx(oxygen_list[i]).GetNeighbors()
        oxygen_atom = big_mol.GetAtomWithIdx(oxygen_list[i])
        neighbor_carbon_idx = [neighbor.GetIdx() for neighbor in oxygen_atom.GetNeighbors()]
    
    # Add neighbor indices to the list #
        list_carbon_idx.append(neighbor_carbon_idx)
        
    # Get nitrogen's index #
    nitrogen_idx = big_mol.GetSubstructMatches(Chem.MolFromSmiles("N"))[0][0]
    
############################# Create the double bond to form a pyrrole pattern #########################################

    ## Get all index of the neighbors of the C1 ##
    C1_neighbor_list = []
    for bond in big_mol.GetBonds():
        if bond.GetBeginAtomIdx() == list_carbon_idx[0][0]:
            C1_neighbor_list.append(bond.GetEndAtomIdx())
        elif bond.GetEndAtomIdx() == list_carbon_idx[0][0]:
            C1_neighbor_list.append(bond.GetBeginAtomIdx())

    ## Get all index of the neighbors of the C2 ##
    C2_neighbor_list = []
    for bond in big_mol.GetBonds():
        if bond.GetBeginAtomIdx() == list_carbon_idx[1][0]:
            C2_neighbor_list.append(bond.GetEndAtomIdx())
        elif bond.GetEndAtomIdx() == list_carbon_idx[1][0]:
            C2_neighbor_list.append(bond.GetBeginAtomIdx())

    ## Get all index of the second neighbor of the C1 ##
    C1_neighbors_neighbor = []
    for i in range(len(C1_neighbor_list)):
        for bond in big_mol.GetBonds():
            if bond.GetBeginAtomIdx() == C1_neighbor_list[i]:
                C1_neighbors_neighbor.append(bond.GetEndAtomIdx())
            elif bond.GetEndAtomIdx() == C1_neighbor_list[i]:
                C1_neighbors_neighbor.append(bond.GetBeginAtomIdx())

    ## Get all index of the second neighbor of the C2 ##
    C2_neighbors_neighbor = []
    for i in range(len(C2_neighbor_list)):
        for bond in big_mol.GetBonds():
            if bond.GetBeginAtomIdx() == C2_neighbor_list[i]:
                C2_neighbors_neighbor.append(bond.GetEndAtomIdx())
            elif bond.GetEndAtomIdx() == C2_neighbor_list[i]:
                C2_neighbors_neighbor.append(bond.GetBeginAtomIdx())  

    ## Find the right neighboring carbon by comparing the index list of the second neighbor of one of the two carbons and the first of the other ##

    C1_neighbors_neighbor_set = set(C1_neighbors_neighbor)
    C1_neighbor_set = set(C1_neighbor_list)
    C2_neighbors_neighbor_set = set(C2_neighbors_neighbor)
    C2_neighbor_set = set(C2_neighbor_list)
    C1_idx = C1_neighbors_neighbor_set.intersection(C2_neighbor_set)
    C2_idx = C2_neighbors_neighbor_set.intersection(C1_neighbor_set)
    
    ## Get the two single bonds index that wil be replace by two double bond ##
    double_bond_idx1 = big_mol.GetBondBetweenAtoms(list_carbon_idx[0][0], list(C2_idx)[0]).GetIdx()
    double_bond_idx2 = big_mol.GetBondBetweenAtoms(list_carbon_idx[1][0], list(C1_idx)[0]).GetIdx()
    
    ## Replace the single bond by double bond ##
    big_mol.GetBondWithIdx(double_bond_idx1).SetBondType(Chem.rdchem.BondType.DOUBLE)
    big_mol.GetBondWithIdx(double_bond_idx2).SetBondType(Chem.rdchem.BondType.DOUBLE)
    
############################# Molecules ... ASSEMBLE! #########################################    

    ## Modify big_mol to become editable ##
    big_mol_editable = Chem.EditableMol(big_mol)
    
    ## the nitrogen of start2 and the neighbor's oxygen carbon ##
    big_mol_editable.AddBond(nitrogen_idx, list_carbon_idx[0][0], order=Chem.rdchem.BondType.SINGLE)
    big_mol_editable.AddBond(nitrogen_idx, list_carbon_idx[1][0], order=Chem.rdchem.BondType.SINGLE)
    
    ## Remove excess oxygen atoms ##
    big_mol_editable.RemoveAtom(oxygen_list[0])
    big_mol_editable.RemoveAtom(oxygen_list[1]-1)
    
###################################### Last Treatment ###################################################
        
    finish_mol1 = big_mol_editable.GetMol()
    Chem.SanitizeMol(finish_mol1)
    finish_mol_smiles = Chem.MolToSmiles(finish_mol1)
    
#########################################################################################

    return finish_mol1


# <span style="font-size:14pt; text-decoration:underline;"> Bisulfite Addition </span>

# In[15]:


def bisulfite_add_compounds1_function(input_smiles):
     ############################# Take the molec apart #########################################

    ## set the variable ##
    start_mol = Chem.MolFromSmiles(input_smiles)
    start_pattern = Chem.MolFromSmiles("C=O")
    finish_pattern = Chem.MolFromSmiles("C(O)(S(=O)(=O)[O-].[Na+])")
    
    ## slicing the molecules to find R1 and R2 ##
    R12 = Chem.DeleteSubstructs(start_mol, start_pattern)
    R12_smiles = Chem.MolToSmiles(R12)

    # slicing for R1 and R2 #
    R_list = R12_smiles.split(".")
    R1_smiles = R_list[0]
    R1 = Chem.MolFromSmiles(R1_smiles)
    R2_smiles = R_list[1]
    R2 = Chem.MolFromSmiles(R2_smiles)

    ## find R1, R2 and start idx ##
    match_start_mol_pattern = start_mol.GetSubstructMatch(start_pattern)
    R1_idx = 0
    R2_idx = len(R1.GetAtoms())
    finish_idx = len(R12.GetAtoms()) + len(finish_pattern.GetAtoms()) - 3

    ############################# Molecules ... ASSEMBLE! #########################################

    ### String Assemble ##
    big_smiles = R12_smiles + "." + Chem.MolToSmiles(finish_pattern)
    big_smiles_mol = Chem.MolFromSmiles(big_smiles)
    big_smiles_mol_editable = Chem.EditableMol(big_smiles_mol)
    #map_big_smiles = mol_with_atom_and_molecule_index(big_smiles_mol)
    
    ## Making bond ##
    big_smiles_mol_editable.AddBond(R1_idx, finish_idx, order=Chem.rdchem.BondType.SINGLE)
    big_smiles_mol_editable.AddBond(R2_idx, finish_idx, order=Chem.rdchem.BondType.SINGLE)
    finish_mol1 = big_smiles_mol_editable.GetMol()
    
    ## Sanitize ##
    Chem.SanitizeMol(finish_mol1)
    
    finish_mol_smiles = Chem.MolToSmiles(finish_mol1)

    return finish_mol1


# <span style="font-size:14pt; text-decoration:underline;"> Cyanohydrine Reduction </span>

# In[16]:


def reduction_cyanohydrine_function(input_smiles):
     ############################# Take the molec apart #########################################

    ## set the variable ##
    start_mol = Chem.MolFromSmiles(input_smiles)
    start_pattern = Chem.MolFromSmiles("C=O")
    finish_pattern = Chem.MolFromSmiles("C(O)(CN)")
    
    ## slicing the molecules to find R1 and R2 ##
    R12 = Chem.DeleteSubstructs(start_mol, start_pattern)
    R12_smiles = Chem.MolToSmiles(R12)

    # slicing for R1 and R2 #
    R_list = R12_smiles.split(".")
    R1_smiles = R_list[0]
    R1 = Chem.MolFromSmiles(R1_smiles)
    R2_smiles = R_list[1]
    R2 = Chem.MolFromSmiles(R2_smiles)

    ## find R1, R2 and start idx ##
    match_start_mol_pattern = start_mol.GetSubstructMatch(start_pattern)
    R1_idx = 0
    R2_idx = len(R1.GetAtoms())
    finish_idx = len(R12.GetAtoms()) + len(finish_pattern.GetAtoms()) - 2

############################# Molecules ... ASSEMBLE! #########################################

    ### String Assemble ##
    big_smiles = R12_smiles + "." + Chem.MolToSmiles(finish_pattern)
    big_smiles_mol = Chem.MolFromSmiles(big_smiles)
    big_smiles_mol_editable = Chem.EditableMol(big_smiles_mol)
    
    ## Making bond ##
    big_smiles_mol_editable.AddBond(R1_idx, finish_idx, order=Chem.rdchem.BondType.SINGLE)
    big_smiles_mol_editable.AddBond(R2_idx, finish_idx, order=Chem.rdchem.BondType.SINGLE)
    finish_mol1 = big_smiles_mol_editable.GetMol()
    
    ## Sanitize ##
    Chem.SanitizeMol(finish_mol1)
    finish_mol_smiles = Chem.MolToSmiles(finish_mol1)
    
    return finish_mol1







