#!/usr/bin/env python
# coding: utf-8

# In[8]:


#get_ipython().run_line_magic('run', './Reaction_Functions_python.py')
#get_ipython().run_line_magic('run', './Dictionnaries_python.py')
from ppc_froi_project.Dictionnaries_python import *
from ppc_froi_project.Reaction_Functions_python import *

# In[9]:


##All the Start patterns that permit to identify the reaction

Start = ["CC(=O)C", "OCCO", "OCCCO" , "NN", "OCCCCC(=O)CCCCO", "CC(O)(O)C", "c1ccc2c(c1)c[nH]c2CCN", 
         "CC(OC)(OC)C", "CN", "C#N", "CC(=O)CCC(=O)C", "CC(=O)", "C1CCC=CC(=O)1", "CC(=O)OC",
         "CC(=O)N(C)C", "C[C@@H](O)CC(=O)C", "c1ccccc1O", "CC(=O)[Si](CC)(CC)CC", "CC=CC(=O)C", 
         "P(c1ccccc1)(c1ccccc1)(c1cccccc1)CC", "P(=O)(OC)(OC)C", "P(OC)(OC)(OC)", "CCI", "CCBr",
         "CCF", "CCCl", "P(OCC(F)(F)F)(OCC(F)(F)F)(=O)CC(=O)OC" ]

Start2 = ["N=O", "NO", "OCCO", "OCCCO", "NC", "CN", "NN", "CC(OC)(OC)C", "c1ccc2c(c1)c[nH]c2CCN", 
          "CC=O", "c1ccccc1O", "CC(=O)C", "CC=CC(=O)C","CC(O)C", "P(c1ccccc1)(c1ccccc1)(c1cccccc1)CC", 
          "P(=O)(OC)(OC)C", "CCI", "CCBr", "CCF", "CCCl", "P(OCC(F)(F)F)(OCC(F)(F)F)(=O)CC(=O)OC"]


##All the Finish patterns tht permit to identify the reactions

Finish = ["CC(O)C", "CC(O)(O)C","CC(O)(OC)C", "CC(OC)(OC)C","C1(C)(C)OCCO1", 
          "C1(C)(C)OCCCO1","CC(=NC)C","CC(=NO)C", "CC(=NN)C", "O1CCCC(C1)1OCCCC1", 
          "C1(C)(C)OC(C)C(C)O1","n1CC2=C(C(C)(C)1)NC3=C2C=CC=C3", "NCCS(=O)(=O)[O-].[Na+]", 
          "CC(O)(S(=O)(=O)[O-].[Na+])C", "CC(N)(C#N)C", "CC(O)(CN)C", "CC(O)(COO)C", "CC(O)(C#N)C",
          "Cn1c(C)ccc1(C)", "CC([O-])=CC", "CC(O)C", "CCO", "C1CCC=CC(O)1", "CCN(C)C", "CC(NC)C", "CCC",
          "CC(=O)OCC", "CC[C@H](O)C[C@H](C)OC(=O)C", "c1ccccc1C(=O)O", "c1ccccc1C(=O)Oc1ccccc1", 
          "CCC(=O)C(O[Si](CC)(CC)CC)C", "CC(=O)C(C)CC(=O)C", "CC=C(C)C", "CC=CC", "P(=O)(OC)(OC)CC",
          "C/C=C\C(=O)OC", "C=O", "c1ccccc1CO", "CO" ]

## All the available reaction conditions 

Condition = [ "Base and Nucleophile", "Acid and Nucleophile", "H2O and acid", " H+ and 1equiv MeOH", 
             "H+ and 2equiv MeOH","Acid", "pH>6", "4<pH<6", "Base", "OCCS(=O)(=O)[O-].[Na+]", 
             "S(=O)(=O)[O-].[Na+]", "NH3, HCN", "HCN, LiAlH4", "H+,H2O", "HCN", "C2OOH", "Base", "H-",
             "NaBH4, -78°C", "NaBH4, CeCl3, MeOH", "LiAlH4, 2 équiv.", "LiAlH4",  "H+, NaBH3CN", "NH2NH2, KOH, 180°C",
             "Zn(Hg), HCl, reflux", "Al(OR)3", "SmI2, THF, -10°C", "KOH, reflux", "CN-, EtOH, H2O",
             "KCN (0,3 equiv), 18-C-6 (0,1 equiv)", "Catalyst: CN", "Al(iPr)3", "KHMDS, 18-Crown-6, -78°C"]

## All the reactions names 

Reaction = [ "wolf_kishner_huang_reduction", "nucleophilic_attack2", "hydration", "Acetalisation1", "Acetalisation2", "Cyclic_Acetal_Formation1",
            "Cyclic_Acetal_Formation2", "Imine_Formation1", "Imine_Formation2", "Oxine_formation", "Hydrazone_formation", 
            "Spiro_Acetal_Formation", "Transacetalisation", 
            "carbonyles_deprotonation", "heterocycle_synthesis", "addition_cyanide", "hydrolysis_cyanohydrine", 
            "reduction_cyanohydrine", "strecker_reaction", "bisulfite_add_compounds1", "bisulfite_add_compounds2", "general_carbonyl_reduction", 
            "aldehyde_reduction", "luche_reaction", "ester_reduction", "amide_reduction", "iminie_reduction", "clemmesen_reduction", 
            "tischenko_reaction", "tischenko_evans_reaction", "canizzaro_reaction", "benzoin_reaction", "cross_benzoin_reaction", "stetter_reaction", 
            "mpv_reaction", "wittig_reaction", "horner_wadsworth_emmons_rxn", "abruzov_reaction_iode", "abruzov_reaction_brome", 
            "abruzov_reaction_fluor", "abruzov_reaction_chlore", "still_grenari_olefination", "pinacol", "nucleophilic_attack1",  "Pictet_Spengler_Reaction"
                ]

## Link between the reaction and its corresponding function

function_reaction = {"wolf_kishner_huang_reduction": WolfKishner_function, "Oxine_formation": OxineFormation_function, 
                     "hydration": Hydration_function, "bisulfite_add_compounds1": bisulfite_add_compounds1_function, 
                     "strecker_reaction": strecker_reaction_function, "reduction_cyanohydrine": reduction_cyanohydrine_function, "pinacol": pinacol_function,
                     "heterocycle_synthesis": heterocycle_synthesis_function, "hydrolysis_cyanohydrine":hydrolysis_cyanohydrine_function,
                     "Imine_Formation2": Imine_Formation2_function, "addition_cyanide": addition_cyanide_function, 
                     "carbonyles_deprotonation": carbonyles_deprotonation_function}


# In[ ]:




