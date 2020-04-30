#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  6 12:12:55 2019

@author: N. Youssef 

Return contact map given pdb structure
"""
import numpy as np
import Bio.PDB
import helpers

pdb_code = '2ppn'

#load contact map # 
with open("../contact_maps/" + protein + "_DICT.txt", "rb") as file:
    ContactMapNS = pickle.load(file)
   
###############################################################################
# energy functions 
###############################################################################
kt = 0.6
def fitness(seq, contact_map):
    Gns = helpers.getG(seq, contact_map) #G in Native structure
    Galt =  helpers.getGalt(seq) #G in alt structures
    dG = Gns - np.mean(Galt) + (np.var(Galt) / (2*kt)) +(kt*np.log(3.4**len(seq))) #delta G
    fitness = np.exp(-1.*dG/kt)/(1+ np.exp(-1.*dG/kt)) #Prob of folding = fitness
    return(fitness, Gns, Galt, dG)
    

def simple_fitness(site, aa, Gns_old, Galt_old, seq, contact_map):
    new_seq = list(seq)
    new_seq[site] = aa
    
    Gns = helpers.simple_getG(site, Gns_old, new_seq, seq, contact_map) #G in Native structure
    Galt =  helpers.simple_getGalt(site, Galt_old, new_seq, seq) #G in alt structures
    
    dG = Gns - np.mean(Galt) + (np.var(Galt) / (2*kt)) +(kt*np.log(3.4**len(seq))) #delta G
    fitness = np.exp(-1.*dG/kt)/(1+ np.exp(-1.*dG/kt)) #Prob of folding = fitness
    return(fitness, np.mean(Galt))
    
###############################################################################
# Functions to evolve 
###############################################################################
def site_specific_fit(seq, contact_map):
	'''
		Calculates the site specific fitness landscape at all sites 
	'''
    ss_fit = np.zeros((len(seq), 20))
    ss_Galt = np.zeros((len(seq), 20))
    fit, Gns_old, Galt_old, dG = fitness(seq, contact_map)
    for site in range(len(seq)):
        for aa in range(20):
            if seq[site] == aa:
                ss_fit[site, aa] = -1
            else:
                new_fit, Gbar = simple_fitness(site, aa, Gns_old, Galt_old, seq, contact_map)
                ss_fit[site, aa] = new_fit
                ss_Galt[site, aa] = Gbar
    return(ss_fit, ss_Galt, fit)


def higher_fitness(ss_fit, seq, current_fit):
    idx = np.random.randint(len(np.where(ss_fit > current_fit)[0]))
    site =  np.where(ss_fit > current_fit)[0][idx]
    aa   =  np.where(ss_fit > current_fit)[1][idx]
    return(site, aa)
    
    
###############################################################################
if pdb_code == '1qhw':
    NS_Seq = 'STLRFVAVGDWGGVPNAPFHTAREMANAKEIARTVQIMGADFIMSLGDNFYFTGVHDANDKRFQETFEDVFSDRALRNIPWYVLAGNHDHLGNVSAQIAYSKISKRWNFPSPYYRLRFKVPRSNITVAIFMLDTVMLCGNSDDFVSQQPEMPRDLGVARTQLSWLKKQLAAAKEDYVLVAGHYPIWSIAEHGPTRCLVKNLRPLLAAYGVTAYLCGHDHNLQYLQDENGVGYVLSGAGNFMDPSVRHQRKVPNGYLRFHYGSEDSLGGFTYVEIGSKEMSITYVEASGKSLFKTSLPRRP'
elif pdb_code == '2ppn':
    NS_Seq = 'GVQVETISPGDGRTFPKRGQTCVVHYTGMLEDGKKFDSSRDRNKPFKFMLGKQEVIRGWEEGVAQMSVGQRAKLTISPDYAYGATGHPGIIPPHATLVFDVELLKLE'
elif pdb_code == '1pek':
    NS_Seq = 'AAQTNAPWGLARISSTSPGTSTYYYDESAGQGSCVYVIDTGIEASHPEFEGRAQMVKTYYYSSRDGNGHGTHCAGTVGSRTYGVVKKTQLFGVKVLDDNGSGQYSTIIAGMDFVASDKNNRNCPKGVVASLSLGGGYSSSVNSAAARLQSSGVMVAVAAGNNNADARNYSPASEPSVCTVGASDRYDRRSSFSNYGSVLDIFGPGTSILSTWIGGSTRSISGTSMATPHVAGLAAYLMTLGKTTAASACRYIADTANKGDLSNIPFGTVNLLAYNNYQA'
    
NS = helpers.AminoAcid_to_AAindex(NS_Seq)
seq = list(NS)
max_fit = -0.01
current_fit = site_specific_fit(seq, NS_map)[2]

while current_fit < 0.99:
    ss_fit, ssGbar, current_fit = site_specific_fit(seq, NS_map)    
    if np.max(ss_fit) < current_fit:        
    # if none of the single step amino acid changes will increase fitness 
    # then randomly choose 20 sites and change them to higher fit aa
        sites = np.random.choice(len(seq), 30)
        new_seq = list(seq)
        for site in sites:
            aa = np.where(ss_fit[site] == ss_fit[site].max())[0][0]
            new_seq[site] = int(aa)
        seq = list(new_seq)        
    
    else:
        site, aa = higher_fitness(ss_fit, seq, current_fit)
        new_seq = list(seq)
        new_seq[site] = int(aa)
        seq = list(new_seq)
            
            
    print(helpers.aa_index_to_aa(seq))
    print(str(current_fit) + '\n\n')
    print(site, aa)
    



