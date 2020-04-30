#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 13:15:13 2019

@author: Noory333

Calculate site-specific fitness landscape for each sequence in each trial 
and calculate Kn and Ln for S-SD model
and write out Average site-specific fitness landscape to use to generate under S-SI
"""
import numpy as np
import helpers
from scipy import linalg
import pickle

protein = '1pek'
path_to_alignments = "../Alignments/"+ protein +"_Ne2_S-SD/"
path_to_results = "../Results/dNdS/"+ protein +"_Ne2_S-SD_Kn_Ln/"
path_to_ssFit = "../Results/site_specific_fitness/"+ protein +"_Ne2_S-SD_ssFit/"

#load contact map #
with open("../contact_maps/" + protein + "_DICT.txt", "rb") as file:
    ContactMapNS = pickle.load(file)


###############################################################################
####          functions to calculate ssfit for extant seqs                 ####
###############################################################################
def ReadSequences(InFile, return_aaSeqs):
    '''
        Reads in an alignment file and returns and array of the sequences (
        return_aaSeqs = 0 => returns CDN idx sequences 
        return_aaSeqs = 1 => returns aa idx sequences
    '''
    with open(InFile) as file:
        line = file.readlines()
    
    Seqs = np.zeros((num_taxa, length))
    for l,s in enumerate(range(2, len(line),2)):
        #print(s,l)
        CodonSeq = line[s]
        Idx = helpers.Codon_to_Cindex(CodonSeq)
        if return_aaSeqs:
            Idx = helpers.Codon_to_AA(Idx)
        Seqs[l] = Idx
    return(Seqs)
    

def ssFitness(SEQ, Site):
    """
        Calculates the site-specific fitness values for all codons given the rest of the sequence is held constant
    """
    Pold, GNSold, GALTold, dGold = helpers.Fitness(SEQ, 1, ContactMapNS)
    F = []
    for codon in range(61):
        SEQnew = list(SEQ); SEQnew[Site] = codon
        Pnew = helpers.simple_Fitness(SEQnew, Site, GNSold, GALTold, SEQ, ContactMapNS)
        F.append(Pnew)
    return(F)
    

def seq_ssFitness(SEQ):
	'''
		Calculates site-specific fitness for all site in the protein
	'''
    ssCDN_F = np.zeros((length, num_taxa, 61))
    for site in range(length):
        for taxa in range(num_taxa):
            seq = CDN_Seqs[taxa]
            ssCDN_F[site, taxa] = ssFitness(seq, site)
    return(ssCDN_F)


###############################################################################
####                    functions to calculate dN/dS                       ####
###############################################################################
def Q_matrix(F, GTR, Neff):
	'''
    creates a 61x61 instantaneous rate matrix based on the fitness vector (F)
    and the mutation model (GTR) and Neff 
	'''
	Q = np.zeros(( 61,61))
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
	return Kn, Ln

def seq_Kn_Ln(ssCDN_F):
    KN = np.zeros((num_taxa,length)); LN = np.zeros((num_taxa,length))
    for site in range(length):
            for taxa in range(num_taxa):
                F = ssCDN_F[site,taxa]
                
                #calculate stationary freq
                pi = [p_neutral[x] * np.exp(2*Neff*F[x]) for x in range(0,61) ]
                pi = pi/np.sum(pi)
                
                #calculate transition matrix 
                Q = Q_matrix(F, GTR, Neff)
                
                #calculate dN/dS
                Kn, Ln = expected_dN_dS(pi, Q, M, IN)
                KN[taxa, site] = Kn
                LN[taxa, site] = Ln
    return(KN, LN)

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

###############################################################################
####                    calculate E[dN/dS] for S-SD                        ####
###############################################################################
for trial in range(1,51):
    print(trial)
    CDN_Seqs = ReadSequences(path_to_alignments + 'seqfile'+ str(trial) + '.txt', 0)
    ssCDN_F = seq_ssFitness(CDN_Seqs)
    pickle.dump(ssCDN_F,open(path_to_ssFit + "ssFit_seqfile"+ str(trial) + ".pkl", 'wb'))

    Kn, Ln = seq_Kn_Ln(ssCDN_F)
    np.savetxt(path_to_results + "Kn_" + str(trial)+ ".csv", Kn, delimiter=" ")
    np.savetxt(path_to_results + "Ln_" + str(trial)+ ".csv", Ln, delimiter=" ")

    
    
##### Write average site-specific fitness to use for S-SI simulations #### 
path_to_ssFit = "../Results/site_specific_fitness/"

for trial in range(1, 51):
    full_ssFit = pickle.load(open(path_to_ssFit + protein +'_Ne2_S-SD_ssFit/ssFit_seqfile'+str(trial)+'.pkl', 'rb'))
    ssFit = np.mean(full_ssFit, axis = 1)
    np.savetxt(path_to_ssFit + protein +"_Ne2_S-SD_AvgssFit/ssFit_seqfile" + str(trial)+ ".csv", ssFit.T, delimiter=" ")

