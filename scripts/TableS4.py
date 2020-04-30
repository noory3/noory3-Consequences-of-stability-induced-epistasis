#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 09:07:38 2020

@author: N. Youssef 

Testing for difference in MLES means across generative models 
"""
import pandas as pd 
import statsmodels.api as sm
from statsmodels.formula.api import ols
import scikit_posthocs as sp
import helpers
import numpy as np

protein = '2ppn'


rst_path = '../Results/Inference/'
path_to_alignments = "../Alignments/"

num_taxa = 14
num_sites = 300 

def taxa_sites(InFile):
    with open(InFile) as file:
        line = file.readlines()

    num_taxa = int(line[0].strip().split(' ')[0])
    num_sites = int(int(line[0].strip().split(' ')[-1])/3)
    return(num_taxa, num_sites)
    
def ReadSequences(InFile):
    '''
        Reads in an alignment file and returns and array of the sequences (AA index)
    '''
    with open(InFile) as file:
        line = file.readlines()
        
    Seqs = np.zeros((num_taxa, num_sites))
    for l,s in enumerate(range(2, len(line),2)):
        CodonSeq = line[s]
        CodonIdx = helpers.Codon_to_Cindex(CodonSeq)
        AAIdx = helpers.Codon_to_AA(CodonIdx)
        Seqs[l] = AAIdx
    return(Seqs)

def CountAAperSite(Seqs):
    '''
        Counts the number of amino acids realized at a site 
    '''
    AACounts = []
    for site in range(num_sites):
        AACounts.append(len(np.unique(Seqs[:,site])))
    return (AACounts)

def Proportions(AACounts):
    '''
        Calculates the proprtion of sites with specific AA counts
    '''
    Prop = np.zeros(num_taxa)
    for i in range(1,num_taxa):
        Prop[i-1] = AACounts.count(i)
    Prop = Prop/num_sites
    return(Prop)

    
def FullSimProp(InFile):
    '''
        Returns full proprtion vector for simulations
    '''
    SimProp = np.zeros((50,num_taxa))
    infile = InFile
    for i in range(1,50+1):
        infile += "/seqfile" + str(i) + ".txt"
        Prop = Proportions(CountAAperSite(ReadSequences(infile)))
        SimProp[i-1] = Prop
        infile = InFile
    return SimProp

def num_aa_per_site(protein):
    global num_taxa
    global num_sites
    num_taxa, num_sites = taxa_sites(path_to_alignments + protein +"_numbered.txt")
        
    if protein == '2ppn':
        num_sites -= 1 
    elif protein == '1pek':
        num_sites = 279
        
    CSI = FullSimProp(path_to_alignments + protein + "_Ne2_C-SI")
    SSI = FullSimProp(path_to_alignments + protein + "_Ne2_S-SI")
    SSD = FullSimProp(path_to_alignments + protein + "_Ne2_S-SD")
    
    return(CSI, SSI, SSD)



#%% Load M3 k=2 data 
CSI      = pd.read_csv(rst_path + protein + '_Ne2_C-SI_M3_k2.csv', sep = " ")
SSI      = pd.read_csv(rst_path + protein + '_Ne2_S-SI_M3_k2.csv', sep = " ")
SSD      = pd.read_csv(rst_path + protein + '_Ne2_S-SD_M3_k2.csv', sep = " ")

# create dataframe
omega_0 = list(CSI.w0) + list(SSI.w0) + list(SSD.w0)
omega_1 = list(CSI.w1) + list(SSI.w1) + list(SSD.w1)

model = (['CSI']*50) + (['SSI']*50)  + (['SSD']*50)

data = pd.DataFrame({'model': model, 'omega_0' : omega_0, 'omega_1': omega_1})

# check if any of the means is significantly different from the rest 
lm = ols('omega_0 ~ model', data = data).fit()
table = sm.stats.anova_lm(lm)
print(table)

# post-hoc test to see if CSI mean is significantly higher 
print('omega_0')
print(sp.posthoc_ttest(data, val_col='omega_0', group_col='model', p_adjust='bonferroni'))

print('\nomega_1')
print(sp.posthoc_ttest(data, val_col='omega_1', group_col='model', p_adjust='bonferroni'))


#%% Load number of amino acids per site
CSI, SSI, SSD = num_aa_per_site(protein)
one = list(CSI[:,0]) + list(SSI[:,0]) + list(SSD[:,0])
five = list(np.sum(CSI[:,4:], axis = 1)) + list(np.sum(SSI[:,4:], axis = 1)) + list(np.sum(SSD[:,4:], axis = 1))

model = (['CSI']*50) + (['SSI']*50)  + (['SSD']*50)
data = pd.DataFrame({'model': model, 'one' : one, 'five': five})


print('# aa per site = 1')
#anova
lm = ols('one ~ model', data = data).fit()
table = sm.stats.anova_lm(lm)
print(table)
#post-hoc
print(sp.posthoc_ttest(data, val_col='one', group_col='model', p_adjust='bonferroni'))

print('\n # aa per site > = 5')
#anova
lm = ols('five ~ model', data = data).fit()
table = sm.stats.anova_lm(lm)
print(table)
#post-hoc
print(sp.posthoc_ttest(data, val_col='five', group_col='model', p_adjust='bonferroni'))