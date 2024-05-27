#!/usr/bin/env python
# coding: utf-8

# In[7]:


#Reactions dictionnaries

## These are the reactions recognized by our code 
## Most of the reactions have a corresponding function that permits to return the final product
## However, for reactions that lack a specific function, the code simply returns the name of the reaction, offering valuable assistance to the user

nucleophilic_attack1 = { "Name":"Nucleophilic attack on basic conditions",
                        "Start": "CC(=O)C", "Start2": "" , "Condition" : "Base and Nucleophile",  "Finish":"CC(O)C", "Finish2":""} 

nucleophilic_attack2 = { "Name": "Nucleophilic attack on acidic conditions",
                        "Start": "CC(=O)C", "Start2": "" , "Condition" : "Acid and Nucleophile",  "Finish":"CC(O)C", "Finish2":""}

hydration = { "Name":"Hydration of Carbonyls" ,
             "Start": "CC(=O)C", "Start2": "" , "Condition" : "H2O and acid",  "Finish":"CC(O)(O)C", "Finish2":""} 

Acetalisation1 = { "Name":"Acetalisation to form hemiacetal ",
             "Start": "CC(=O)C", "Start2": "" , "Condition" : " H+ and 1equiv MeOH",  "Finish":"CC(O)(OC)C", "Finish2":""} 

Acetalisation2 = { "Name":"Acetalisation to form acetal",
             "Start": "CC(=O)C", "Start2": "" ,"Condition" : " H+ and 2equiv MeOH ",  "Finish":"CC(OC)(OC)C", "Finish2":""} 

Cyclic_Acetal_Formation1= { "Name":"Formation of cyclic acetal from diols addition with 2 carbons",
                           "Start": "CC(=O)C" , "Start2":"OCCO", "Condition" : "Acid",  "Finish":"C1(C)(C)OCCO1", "Finish2":""} 

Cyclic_Acetal_Formation2= { "Name":"Formation of cyclic acetal from diols addition with 3 carbons ",
                           "Start": "CC(=O)C" , "Start2":"OCCCO", "Condition" : "Acid",  "Finish":"C1(C)(C)OCCCO1", "Finish2":""} 

Imine_Formation1 = { "Name": "Addition of amines on basic conditions" ,
                    "Start": "CC(=O)C" , "Start2":"NC","Condition" : "pH>6",  "Finish":"CC(=NC)C", "Finish2":""}

Imine_Formation2 = { "Name": " Addition of amines on acidic conditions",
                    "Start": "CC(=O)C" ,"Start2":"NC", "Condition" : "4<pH<6",  "Finish":"CC(=NC)C" , "Finish2":""}

Oxine_formation = { "Name": "Addition of nucleophilic nitrogen to form an oxine" ,
                  "Start": "CC(=O)C" , "Start2":"NO", "Condition" : "",  "Finish":"CC(=NO)C" , "Finish2":""}

Hydrazone_formation = { "Name": "Addition of nucleophilic nitrogen to form a hydrozone" ,
                      "Start": "CC(=O)C" , "Start2":"NN", "Condition" : "",  "Finish":"CC(=NN)C" , "Finish2":""}

Spiro_Acetal_Formation = { "Name":"Formation of Spiroacetal with 2 cycles",
                          "Start": "OCCCCC(=O)CCCCO", "Start2": "" , "Condition" : "Acid",  "Finish":"O1CCCC(C1)1OCCCC1" , "Finish2":""}

Transacetalisation = { "Name": "Transacetalisation to obtain a spiro"  ,
                      "Start": "CC(O)(O)C" , "Start2":"CC(OC)(OC)C", "Condition" : "Acid",  "Finish":"C1(C)(C)OC(C)C(C)O1", "Finish2":"" }

Pictet_Spengler_Reaction = { "Name": "Pictet Spengler Reaction"  ,
                            "Start": "CC(=O)C" ,"Start2":"c1ccc2c(c1)c[nH]c2CCN", "Condition" : "Base",
                            "Finish":"n1CC2=C(C(C)(C)1)NC3=C2C=CC=C3" , "Finish2":""}

bisulfite_add_compounds2 = {"Finish": "NCCS(=O)(=O)[O-].[Na+]","Condition" : "OCCS(=O)(=O)[O-].[Na+]", 
                            "Start" : "CN", "Start2": "", "Name": "Bisulfite addition compound Case 2", "Finish2":""}

bisulfite_add_compounds1 = {"Finish": "CC(O)(S(=O)(=O)[O-].[Na+])C", "Condition" : "S(=O)(=O)[O-].[Na+]", 
                            "Start" : "CC(=O)C", "Start2": "", "Name": "Bisulfite addition compound Case 1", "Finish2":""}

strecker_reaction = {"Finish": "CC(N)(C#N)C", "Condition" : "NH3, HCN",
                     "Start" : "CC(=O)C", "Start2": "", "Name": "Strecker Reaction", "Finish2":""}

reduction_cyanohydrine = {"Finish": "CC(O)(CN)C", "Condition" : "HCN, LiAlH4", 
                          "Start" : "CC(=O)C", "Start2": "", "Name": "Reduction of cyanohydrine", "Finish2":""}

hydrolysis_cyanohydrine = {"Finish": "CC(O)(COO)C", "Condition" : "H+,H2O", 
                           "Start" : "CC(=O)C", "Start2": "", "Name": "Hydrolysis of cyanohydrine", "Finish2":""}

addition_cyanide = {"Finish": "CC(O)(C#N)C", "Condition" : "HCN", 
                    "Start" : "N#CC(O)C", "Start2": "", "Name": "Addition of Cyanide", "Finish2":""}

heterocycle_synthesis = {"Finish": "Cn1c(C)ccc1(C)", "Condition" : "C2OOH", 
                         "Start" : "CC(=O)CCC(=O)C", "Start2" : "CN", "Name": "Heterocycle Synthesis", "Finish2":""}

carbonyles_deprotonation = {"Finish": "CC([O-])=CC", "Condition" : "Base", 
                            "Start" : "CC(=O)CC", "Start2": "", "Name": "Deprotonation of Carbonyl", "Finish2":""}

general_carbonyl_reduction = { "Finish" : "CC(O)C", "Condition" : "H-", 
                              "Start" : "CC(=O)C", "Start2": "", "Name": "Reduction of carbonyl", "Finish2":""}

aldehyde_reduction = { "Finish" : "CCO", "Condition" : "NaBH4, -78°C", 
                      "Start" : "CC(=O)", "Start2": "", "Name" : "Reduction of Aldehyde" , "Finish2":""}

luche_reaction = { "Finish" : "C1CCC=CC(O)1", "Condition" : "NaBH4, CeCl3, MeOH", 
                  "Start" : "C1CCC=CC(=O)1", "Start2": "", "Name" : "Luche reaction" , "Finish2":""}

ester_reduction = { "Finish" : "CCO", "Finish2" : "CO", "Condition" : "LiAlH4, 2 équiv.", 
                   "Start" : "CC(=O)OC", "Start2": "", "Name" : "Reduction of ester" , "Finish2":""}

amide_reduction = { "Finish" : "CCN(C)C", "Condition" : "LiAlH4", 
                   "Start" : "CC(=O)N(C)C", "Start2": "", "Name" : "Reduction of amide" , "Finish2":""}

iminie_reduction = { "Finish" : "CC(NC)C", "Condition" : "H+, NaBH3CN", 
                    "Start" : "CC(=O)C", "Start2" : "CN", "Name" : "Reduction of Imine" , "Finish2":""}

wolf_kishner_huang_reduction = { "Finish" : "CCC", "Condition" : "NH2NH2, KOH, 180°C", 
                                "Start" : "CC(=O)C", "Start2": "", "Name" : "Wolf-Kishner-Huang Reduction", "Finish2":"" }

clemmesen_reduction = { "Finish" : "CCC", "Condition" : "Zn(Hg), HCl, reflux",
                       "Start" : "CC(=O)C", "Start2": "", "Name" : "Clemmensen Reduction" , "Finish2":""}

tischenko_reaction = { "Finish" : "CC(=O)OCC", "Condition" : "Al(OR)3", 
                      "Start" : "CC=O", "Start2" : "CC=O", "Name" : "Tischenko Reaction", "Finish2":"" }

tischenko_evans_reaction = { "Finish" : "CC[C@H](O)C[C@H](C)OC(=O)C", "Condition" : "SmI2, THF, -10°C", 
                            "Start" : "C[C@@H](O)CC(=O)C", "Start2" : "CC(=O)", "Name" : "Tischenko-Evans Reaction" , "Finish2":""}

canizzaro_reaction = { "Finish" : "c1ccccc1C(=O)O", "Finish2" : "c1ccccc1CO", "Condition" : "KOH, reflux", 
                      "Start" : "c1ccccc1O", "Start2" : "c1ccccc1O", "Name" : "Cannizzaro Reaction" }

benzoin_reaction = { "Finish" : "c1ccccc1C(=O)Oc1ccccc1", "Condition" : "CN-, EtOH, H2O",
                    "Start" : "c1ccccc1CO", "Start2": "", "Name" : "Benzoin Reaction" , "Finish2":""}

cross_benzoin_reaction = { "Finish" : "CCC(=O)C(O[Si](CC)(CC)CC)C", "Condition" : "KCN (0,3 equiv), 18-C-6 (0,1 equiv)", 
                          "Start" : "CC(=O)[Si](CC)(CC)CC", "Start2" : "CC(=O)", "Name" : "Cross Benzoin Reaction", "Finish2":"" }

stetter_reaction = { "Finish" : "CC(=O)C(C)CC(=O)C", "Condition" : "Catalyst: CN", 
                    "Start" : "CC(=O)", "Start2" : "CC=CC(=O)C", "Name" : "Stetter Reaction", "Finish2":"" }

mpv_reaction = { "Finish" : "CC(O)C", "Finish2" : "C=O", "Condition" : "Al(iPr)3", 
                "Start" : "CC(=O)C", "Start2" : "CC(O)C", "Name" : "Meerwein-Ponndorf-Verley Reaction", "Finish2":"" }

wittig_reaction = { "Finish" : "CC=C(C)C",  "Condition" : "", 
                   "Start" : "CC(=O)C", "Start2" : "P(c1ccccc1)(c1ccccc1)(c1cccccc1)CC", "Name" : "Wittig Reaction" , "Finish2":""}

horner_wadsworth_emmons_rxn = { "Finish" : "CC=CC", "Condition" : "Base", 
                               "Start" : "CC(=O)", "Start2" : "P(=O)(OC)(OC)C", "Name" : "Horner-Wadsworth-Emmons Reaction" , "Finish2":""}

abruzov_reaction_iode = { "Finish" : "P(=O)(OC)(OC)CC", "Condition" : "",
                         "Start" : "P(OC)(OC)(OC)", "Start2" : "CCI",  "Name" : "Abruzov Reaction", "Finish2":"" }

abruzov_reaction_brome = { "Finish" : "P(=O)(OC)(OC)CC", "Condition" : "", 
                          "Start" : "P(OC)(OC)(OC)", "Start2" : "CCBr",  "Name" : "Abruzov Reaction", "Finish2":"" }

abruzov_reaction_fluor = { "Finish" : "P(=O)(OC)(OC)CC", "Condition" : "", 
                          "Start" : "P(OC)(OC)(OC)", "Start2" : "CCF",  "Name" : "Abruzov Reaction" , "Finish2":""}

abruzov_reaction_chlore = { "Finish" : "P(=O)(OC)(OC)CC", "Condition" : "", 
                           "Start" : "P(OC)(OC)(OC)", "Start2" : "CCCl",  "Name" : "Abruzov Reaction", "Finish2":"" }

still_grenari_olefination = { "Finish" : "C/C=C\C(=O)OC", "Condition" : "KHMDS, 18-Crown-6, -78°C", 
                             "Start" : "CC(=O)", "Start2" : "P(OCC(F)(F)F)(OCC(F)(F)F)(=O)CC(=O)OC", "Name" : "Still-Genari Olefination", "Finish2":"" }

pinacol = {"Start": "C(O)C(O)","Start2": "", "Condition" : "H+", "Finish" : "C=O", "Finish2" :"", "Name" : "Pinacol reaction"}
