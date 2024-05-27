#!/usr/bin/env python
# coding: utf-8

# <div style="font-size:24px;"><b>Fonctions & Réactions Organiques I - Practical Programming in Chemistry</b></div>
# 

# <div style="font-size:11pt;"> 
# This project is a digital assistant designed to help students with "Fonctions et Réactions Organiques I" exams. The core function of this code is to input starting materials and conditions, and the tool uses a pre-programmed database to predict and output the resultant chemical compounds</div>
# 

# <span style="font-size:14pt; text-decoration:underline;"> Necessary Imports </span>

# In[1]:

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


# <span style="font-size:14pt; text-decoration:underline;"> Part I: Database </span>
# 
# 

# <div style="font-size:11pt;"> This data section provides an overview of the available chemical reactions, detailing the recognized chemical patterns and the specific conditions under which these reactions occur. This database enables the tool to accurately match user inputs with the appropriate reactions. </div>

# In[2]:


# <span style="font-size:14pt; text-decoration:underline;"> Part II: Reaction Functions </span>

# <div style="font-size:11pt;"> This reaction functions section provides a detailed overview of available reactions and they were developped in order to simulate the correct final compound. These functions enable the tool to accurately react the inital molecule to yield the corresponding final molecule. </div>

# In[3]:


#get_ipython().run_line_magic('run', './Reaction_Functions_python.py')


# <span style="font-size:14pt; text-decoration:underline;"> Part III: Exercice Resolution </span>

# <div style="font-size:11pt;"> _____ Once the tool accurately pair user inputs with corresponding reactions after entering a start molecule and its conditions, the database compares this to its records, identifies the appropriate reaction, and simulates the outcome. This process results in the prediction and visualization of the product molecule, facilitating effective chemical analysis and experimentation. </div>

# In[1]:

from ppc_froi_project.Reaction_Functions_python import *
from ppc_froi_project.Lists_python import *
from ppc_froi_project.Dictionnaries_python import *

## complementary functions
# Function to check if the entered molecule is valid
def validate_smiles(smiles):
    if not smiles:
        return None
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is not None:
        return molecule
    else:
        print("Invalid SMILES string entered.")
        return None

# Function to check for subgroups in start molecules
def subgroup_in_start(mol, Start):
    matched_subgroup = []
    for substr_smiles in Start:
        substructure = Chem.MolFromSmarts(substr_smiles)
        if mol.HasSubstructMatch(substructure):
            matched_subgroup.append(substr_smiles)
            matched_subgroup is not None
    if len(matched_subgroup) == 0:
        print("No Match Found")
        return None 
    else: 
        return matched_subgroup

def subgroup_in_finish(finish_mol, Finish):
    matched_subgroup_f = []
    for substr_smiles in Finish:
        substructure = Chem.MolFromSmarts(substr_smiles)
        if finish_mol.HasSubstructMatch(substructure):
            matched_subgroup_f.append(substr_smiles)
            print(f"Match found for: {substr_smiles}")

    if len(matched_subgroup_f) == 0:
        print("No Match Found")
        return None 
    else:
        return matched_subgroup_f


# Function to check if condition is in the predefined list
def check_for_cond(condition, Condition):
    for i in range(len(Condition)):
        if condition == Condition[i]:
            return True

def check_same_dict_in_list(start, condition, Reaction):
        globals_dict = globals()  # Get the current global symbol table
        matched_subgroup = subgroup_in_start(start, Start)
        for name in Reaction:
            d = globals_dict.get(name)
            #print(f"Checking dictionary '{name}': {d}")
            if d is not None:
                #print(f"Start in dictionary: {matched_subgroup in d['Start']}")
                #print(f"Condition in dictionary: {condition in d['Condition']}"
                for i in range(len(matched_subgroup)):
                      if matched_subgroup[i] in d["Start"] and condition in d["Condition"]:
                        return name  # Return the name of the dictionary
        
        return None  # Return None if the dictionaries are not found
    
def check_same_dict_in_list1(start, start2, Reaction):
    globals_dict = globals()  # Get the current global symbol table
    matched_subgroup = subgroup_in_start(start, Start)
    matched_subgroup2 = subgroup_in_start(start2, Start2)
    for name in Reaction:
        d = globals_dict.get(name)
        #print(f"Checking dictionary '{name}': {d}")
        if d is not None:
            for i in range(len(matched_subgroup)):
                for j in range(len(matched_subgroup2)):
                    if matched_subgroup[i] in d["Start"] and matched_subgroup2[j] in d["Start2"]:  # Update this line to use the correct keys
                        return name  # Return the name of the dictionary
    
    return None  # Return None if the dictionaries are not found

def check_same_dict_in_list3(finish, condition, Reaction):
    globals_dict = globals()  # Access the global scope
    if not matched_subgroup:
        print("No matching subgroup found.")
        return None
    for name in Reaction:
        d = globals_dict.get(name)
        if d is not None:
            for i in range(len(matched_subgroup)):
                if condition in d["Condition"] and matched_subgroup[i] in d["Finish"]:
                    return name
                else:
                    pass

    return None
    
def check_same_dict_in_list4(condition, finish, finish2, Reaction):
    globals_dict = globals()  # Get the current global symbol table
    matched_subgroup = subgroup_in_finish(finish, Finish)
    matched_subgroup2 = subgroup_in_start(finish2, Finish)
    for name in Reaction:
        d = globals_dict.get(name)
        #print(f"Checking dictionary '{name}': {d}")
        if d is not None:
            for i in range(len(subgroup_in_finish)):
                for j in range(len(subgroup_in_start)):
                    if condition in d["Condition"] and matched_subgroup[i] in d["Finish"] and matched_subgroup2[j] in d["Finish2"]:  # Update this line to use the correct keys
                        return name  # Return the name of the dictionary
    
    return None  # Return None if the dictionaries are not found

def check_same_dict_in_list5(start,start2,condition,Reaction):
    globals_dict = globals()  # Get the current global symbol table
    matched_subgroup = subgroup_in_start(start, Start)
    matched_subgroup2 = subgroup_in_start(start2, Start2)
    for name in Reaction:
        d = globals_dict.get(name)
        #print(f"Checking dictionary '{name}': {d}")
        if d is not None:
            for i in range(len(matched_subgroup)):
                for j in range(len(matched_subgroup2)):
                    if matched_subgroup[i] in d["Start"] and matched_subgroup2[j] in d["Start2"] and condition in d["Condition"]:  # Update this line to use the correct keys
                        return name  # Return the name of the dictionary
                        
def apply_function(dictionary_name, dictionary, start_brut):
    func = function_reaction.get(dictionary_name)
    if func:
        start_input_smiles = start_brut # Get the SMILES string from the dictionary
        result = func(start_input_smiles)  # Call the function
        return {'result': result, 'dictionary_name': dictionary_name}  # Return a dictionary with both the result and dictionary name
    else:
        raise ValueError(f"No function defined for dictionary '{dictionary_name}'.")

def apply_function2(dictionary_name, dictionary, start_brut, start2_brut):
    func = function_reaction.get(dictionary_name)
    if func:
        start_input_smiles = start_brut  # Get the SMILES string from the dictionary
        start_input_smiles2 = start2_brut
        result = func(start_input_smiles, start_input_smiles2)  # Call the function
        return {'result': result, 'dictionary_name': dictionary_name}  # Return a dictionary with both the result and dictionary name
    else:
        raise ValueError(f"No function defined for dictionary '{dictionary_name}'.")


# In[2]:


def main():
#for start molecule
    start_brut = input("Initial molecule? Write in SMILES. Else, press Enter")
    start = validate_smiles(start_brut)
    if start:
        if subgroup_in_start(start, Start):
            pass
        else:
            print("The initial molecule you entered cannot be found in our start database.")
#for start molecule 2
    start2_brut = input("Is there a second initial molecule? If not, press Enter")
    start2 = validate_smiles(start2_brut)
    if start2:
        if subgroup_in_start(start2, Start2):  # Pass the molecule object, not the SMILES string
            pass
        else:
            print("The second initial molecule you entered cannot be found in our start database.")
#for conditions
    while True:  # Creating an infinite loop that will continue until explicitly broken out of
        if_condition = input("Conditions? y/n")
        if if_condition.lower() == "y":
            print("Copy paste conditions from the list below")
            print(Condition)
            condition = input("Conditions:")
            if check_for_cond(condition, Condition):
                break  # Exit the loop since the condition is valid
            else:
                print("Condition entered is not in the list. Please try again.")
        elif if_condition.lower() == "n":
            condition = None
            break  # Exit the loop since no conditions are applied
        else:
            print("Invalid input, please enter 'y' for yes or 'n' for no.")
    

    # Handle the final molecule
    finish_brut = input("Final molecule? Write in SMILES. Else, press Enter")
    finish = validate_smiles(finish_brut)
    if finish:
        if subgroup_in_finish(finish, Finish):  # Pass the molecule object, not the SMILES string
            pass
        else:
            print("The final molecule you entered cannot be found in our start database.")
            
    # Handle the second final molecule
    finish2_brut = input("Is there a second final molecule? If not, press Enter")
    finish2 = validate_smiles(finish2_brut)
    if finish2:
        if subgroup_in_finish(finish2, Finish):  # Pass the molecule object, not the SMILES string
            pass
        else:
            print("The final molecule you entered cannot be found in our start database.")

    
    if start and not start2 and condition:
        dictionary_name = check_same_dict_in_list(start, condition, Reaction)
        if dictionary_name: # Check if a dictionary name is found
            result_dict = apply_function(dictionary_name, globals()[dictionary_name], start_brut)
            result = result_dict['result']
            dictionary_name = result_dict['dictionary_name']
            print("Result:", result)
            print("Dictionary name:", dictionary_name)
            print("initial molecule:")
            display(Draw.MolToImage(start))
            print("final molecule:")
            display(Draw.MolToImage(result))
        else:
            print("No matching dictionary found.")
    
    elif start and start2 and condition:
        dictionary_name = check_same_dict_in_list5(start, start2, condition, Reaction)
        if dictionary_name: # Check if a dictionary name is found
            result_dict = apply_function2(dictionary_name, globals()[dictionary_name],start_brut, start2_brut)
            result = result_dict['result']
            dictionary_name = result_dict['dictionary_name']
            print("Result:", result)
            print("Dictionary name:", dictionary_name)
            print("initial molecule:")
            display(Draw.MolToImage(start))
            print("final molecule:")
            display(Draw.MolToImage(result))
        else:
            print("No matching dictionary found.")        
    
    elif start and start2 and not condition:
        dictionary_name = check_same_dict_in_list1(start, start2, Reaction)
        if dictionary_name: # Check if a dictionary name is found
            result_dict = apply_function2(dictionary_name, globals()[dictionary_name],start_brut, start2_brut)
            result = result_dict['result']
            dictionary_name = result_dict['dictionary_name']
            print("Result:", result)
            print("Dictionary name:", dictionary_name)
            print("initial molecule:")
            display(Draw.MolToImage(start))
            print("final molecule:")
            display(Draw.MolToImage(result))
        else:
            print("No matching dictionary found.")

    elif condition and finish:
        print("etape3", type(finish), finish)
        dictionary_name = check_same_dict_in_list3(finish, condition, Reaction)
        print("Dictionary name:", dictionary_name)


    elif condition and finish and finish2:
        dictionary_name = check_same_dict_in_list4(condition, finish, finish2, Reaction)
        print("Dictionary name:", dictionary_name)

if __name__ == "__main__":
    main()