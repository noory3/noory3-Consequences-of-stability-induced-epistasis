#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 10:16:41 2019

@author: Noory333
"""
import numpy as np
import helpers
from scipy import linalg
import pickle

path_to_results = "../Results/"

#%%
protein = '1pek'
path_to_dNdS = path_to_results + "dNdS/" + protein + "_Ne2_C-SI" 
path_to_ssFit = path_to_results + "site_specific_fitness/"+ protein +"_Ne2_C-SI/"

def Q_matrix(F, GTR, Neff):
	'''
    creates a 61x61 instantaneous rate matrix based on the fitness vector (F)
    and the mutation model (GTR) and Neff 
	'''
	Q = np.zeros((61,61))
	for codon1 in range(0,61):
		for codon2 in range(0,61):
			nuc_diff = helpers.NucDiff(helpers.Codon[codon1], helpers.Codon[codon2])
			if len(nuc_diff) == 2:
				Sij = 2*Neff*(F[codon2] - F[codon1])
				n1 = helpers.Nucleotide[nuc_diff[0]]
				n2 = helpers.Nucleotide[nuc_diff[1]]
				if abs(Sij) <= 1e-8:
				    Q[codon1, codon2] = GTR[n1,n2]
				else:   
				    Q[codon1, codon2] = GTR[n1,n2] * (Sij)/(1 - np.exp(-1.*Sij))
	for i in range(0, 61):
		Q[i,i] = -np.sum(Q[i])
	return Q


def IndicatorMatrix():
    '''
    creates indicator matrices specifying if a substitution is synonymous (IS)
    or nonsynonymous (IN)
    '''
    IN = np.zeros((61, 61))
    IS = np.zeros((61, 61))
    for r in range(61):
        for c in range(61):
            diff = helpers.NucDiff(helpers.Codon[r], helpers.Codon[c])
            if len(diff) == 2:
                oldAA = helpers.Codon_AA[helpers.Codon[r]]
                newAA = helpers.Codon_AA[helpers.Codon[c]]
                if oldAA == newAA:
                    IS[r, c] = 1
                else:
                    IN[r, c] = 1
    return (IN, IS)

def expected_dN_dS(pi, Q, M, IN):
	Kn = np.sum(np.multiply(np.dot(np.diag(pi),Q),IN))
	Ln = np.sum(np.multiply(np.dot(np.diag(pi),M),IN))
	dN_dS = Kn / Ln
	return dN_dS, Kn, Ln



if protein == '1qhw':
    pi_nuc = [0.19675, 0.31761, 0.28032, 0.20532] 
    GTR = helpers.MutationMatrix(pi_nuc, 4.49765,1,1,1,1,4.49765)
    num_taxa = int(14)
    length = int(300)

elif protein == '2ppn':
    pi_nuc = [0.19246, 0.24559, 0.29365, 0.26830]
    GTR = helpers.MutationMatrix(pi_nuc, 2.50275 ,1,1,1,1,2.50275)
    num_taxa = int(14)
    length = int(107)

elif protein == '1pek':
    pi_nuc = [0.20853, 0.34561, 0.25835, 0.18750]
    GTR = helpers.MutationMatrix(pi_nuc, 0.90382 ,1,1,1,1, 0.90382)
    num_taxa = int(12)
    length = int(279)
    
    
Neff = int(1e2)
IN , IS = IndicatorMatrix()
M = Q_matrix(np.ones(61), GTR, Neff)

 
#Calculate neutral frequencies based on GTR model 
F_neutral = np.ones((61))
Q_neutral = Q_matrix(F_neutral, GTR, Neff)
P_neutral = linalg.expm( np.multiply(Q_neutral, 40 ) )
p_neutral = P_neutral[0]



#%%############################################################################
####                    calculate E[dN/dS] for  C-SI                       ####
###############################################################################
full_Kn = np.zeros((50,length)); full_Ln = np.zeros((50,length))
for trial in range(1, 51):
    KN = []; LN = []
    print(trial)
    ssFit = np.genfromtxt(path_to_ssFit + 'ssFit_seqfile'+ str(trial) +'.csv', delimiter= ' ')
    for site in range(length):
        F = ssFit[:,site]
        
        #calculate stationary freq
        pi = [p_neutral[x] * np.exp(2*Neff*F[x]) for x in range(0,61) ]
        pi = pi/np.sum(pi)
        
        #calculate transition matrix 
        Q = Q_matrix(F, GTR, Neff)
        
        #calculate dN/dS
        dnds, Kn, Ln = expected_dN_dS(pi, Q, M, IN)
        KN.append(Kn); LN.append(Ln)
    full_Kn[trial-1] = KN
    full_Ln[trial-1] = LN

np.savetxt(path_to_dNdS + "_Kn.csv", full_Kn, delimiter=" ")
np.savetxt(path_to_dNdS + "_Ln.csv", full_Ln, delimiter=" ")
