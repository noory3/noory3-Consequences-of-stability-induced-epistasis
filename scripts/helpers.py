#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 16 15:12:51 2018

@author: N. Youssef

updated helpers functions 
"""

import numpy as np 
import numba as nb
from joblib import Parallel, delayed
import pickle
import csv

location = "../scripts/contact_maps/"
#location = "/Users/nooryoussef/Desktop/Stokes_AntiStokes/Code/Generating/contact_maps/"

############# Contact Maps #############
with open(location + "ALT_DICT.txt", "rb") as file:
    ContactMaps = pickle.load(file)
    
### AMINO ACIDS TO NUMBER INDEX ###
AminoAcid = {'A': 0, 'R': 1, 'N': 2, 'D': 3, 'C': 4, 'Q':5,
            'E': 6, 'G': 7, 'H': 8, 'I': 9, 'L': 10, 'K': 11,
            'M': 12, 'F': 13, 'P': 14, 'S': 15, 'T': 16, 'W': 17, 'Y': 18, 'V': 19}

### NUMBER INDEX TO CODONS ###
Codon = { 0: 'CGT', 1: 'CGC', 2: 'CGA', 3: 'CGG', 4: 'AGA', 5: 'AGG',
              6: 'AAA', 7: 'AAG', 8: 'AAT', 9: 'AAC', 10: 'GAT', 11: 'GAC',
             12: 'GAA', 13: 'GAG', 14: 'CAA', 15: 'CAG', 16: 'CAT',
             17: 'CAC', 18: 'CCT', 19: 'CCC', 20: 'CCA', 21: 'CCG',
             22: 'TAT', 23: 'TAC', 24: 'TGG', 25: 'TCT', 26: 'TCC',
             27: 'TCA', 28: 'TCG', 29: 'AGT', 30: 'AGC', 31: 'ACT',
             32: 'ACC', 33: 'ACA', 34: 'ACG', 35: 'GGT', 36: 'GGC',
             37: 'GGA', 38: 'GGG', 39: 'GCT', 40: 'GCC', 41: 'GCA',
             42: 'GCG', 43: 'ATG', 44: 'TGT', 45: 'TGC', 46: 'TTT',
             47: 'TTC', 48: 'TTA', 49: 'TTG', 50: 'CTT', 51: 'CTC',
             52: 'CTA', 53: 'CTG', 54: 'GTT', 55: 'GTC', 56: 'GTA',
             57: 'GTG', 58: 'ATT', 59: 'ATC', 60: 'ATA'
             }


### CODONS TO AMINO ACIDS ###
Codon_AA =  { "TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
             "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
             "TAT":"Y", "TAC":"Y", 
             "TGT":"C", "TGC":"C", "TGG":"W",
             "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
             "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
             "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
             "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
             "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
             "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
             "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
             "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
             "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
             "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
             "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
             "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G" }

### AA index to synonymous codon index ### 
AA_SYNCDNS =  {0: [39, 40, 41, 42],
 1: [0, 1, 2, 3, 4, 5],
 2: [8, 9],
 3: [10, 11],
 4: [44, 45],
 5: [14, 15],
 6: [12, 13],
 7: [35, 36, 37, 38],
 8: [16, 17],
 9: [58, 59, 60],
 10: [48, 49, 50, 51, 52, 53],
 11: [6, 7],
 12: [43],
 13: [46, 47],
 14: [18, 19, 20, 21],
 15: [25, 26, 27, 28, 29, 30],
 16: [31, 32, 33, 34],
 17: [24],
 18: [22, 23],
 19: [54, 55, 56, 57]}
 
### NUCLEOTIDE TO NUMBER INDEX ###
Nucleotide = {"T":0, "C":1, "G":2, "A":3}


### Contact Poteintials from Miyazawa and Jernigan 1996 ###
MJ85 = np.array(
[[-0.13,  0.43,  0.28,  0.12,  0.00,  0.08,  0.26, -0.07,  0.34, -0.22, -0.01,  0.14,  0.25,  0.03,  0.10, -0.06, -0.09, -0.09,  0.09, -0.10,  0.00],
 [ 0.43,  0.11, -0.14, -0.72,  0.24, -0.52, -0.74, -0.04, -0.12,  0.42,  0.35,  0.75,  0.31,  0.41, -0.38,  0.17, -0.35, -0.16, -0.25,  0.30,  0.00],
 [ 0.28, -0.14, -0.53, -0.30,  0.13, -0.25, -0.32, -0.14, -0.24,  0.53,  0.30, -0.33,  0.08,  0.18, -0.18, -0.14, -0.11,  0.06, -0.20,  0.50,  0.00],
 [ 0.12, -0.72, -0.30,  0.04,  0.03, -0.17, -0.15, -0.22, -0.39,  0.59,  0.67, -0.76,  0.65,  0.39,  0.04, -0.31, -0.29,  0.24,  0.00,  0.58,  0.00],
 [ 0.00,  0.24,  0.13,  0.03, -1.06,  0.05,  0.69, -0.08, -0.19,  0.16, -0.08,  0.71,  0.19, -0.23,  0.00, -0.02,  0.19,  0.08,  0.04,  0.06,  0.00],
 [ 0.08, -0.52, -0.25, -0.17,  0.05,  0.29, -0.17, -0.06, -0.02,  0.36,  0.26, -0.38,  0.46,  0.49, -0.42, -0.14, -0.14,  0.08, -0.20,  0.24,  0.00],
 [ 0.26, -0.74, -0.32, -0.15,  0.69, -0.17, -0.03,  0.25, -0.45,  0.35,  0.43, -0.97,  0.44,  0.27, -0.10, -0.26,  0.00,  0.29, -0.10,  0.34,  0.00],
 [-0.07, -0.04, -0.14, -0.22, -0.08, -0.06,  0.25, -0.38,  0.20,  0.25,  0.23,  0.11,  0.19,  0.38, -0.11, -0.16, -0.26,  0.18,  0.14,  0.16,  0.00],
 [ 0.34, -0.12, -0.24, -0.39, -0.19, -0.02, -0.45,  0.20, -0.29,  0.49,  0.16,  0.22,  0.99, -0.16, -0.21, -0.05, -0.19, -0.12, -0.34,  0.19,  0.00],
 [-0.22,  0.42,  0.53,  0.59,  0.16,  0.36,  0.35,  0.25,  0.49, -0.22, -0.41,  0.36, -0.28, -0.19,  0.25,  0.21,  0.14,  0.02,  0.11, -0.25,  0.00],
 [-0.01,  0.35,  0.30,  0.67, -0.08,  0.26,  0.43,  0.23,  0.16, -0.41, -0.27,  0.19, -0.20, -0.30,  0.42,  0.25,  0.20, -0.09,  0.24, -0.29,  0.00],
 [ 0.14,  0.75, -0.33, -0.76,  0.71, -0.38, -0.97,  0.11,  0.22,  0.36,  0.19,  0.25,  0.00,  0.44,  0.11, -0.13, -0.09,  0.22, -0.21,  0.44,  0.00],
 [ 0.25,  0.31,  0.08,  0.65,  0.19,  0.46,  0.44,  0.19,  0.99, -0.28, -0.20,  0.00,  0.04, -0.42, -0.34,  0.14,  0.19, -0.67, -0.13, -0.14,  0.00],
 [ 0.03,  0.41,  0.18,  0.39, -0.23,  0.49,  0.27,  0.38, -0.16, -0.19, -0.30,  0.44, -0.42, -0.44,  0.20,  0.29,  0.31, -0.16,  0.00, -0.22,  0.00],
 [ 0.10, -0.38, -0.18,  0.04,  0.00, -0.42, -0.10, -0.11, -0.21,  0.25,  0.42,  0.11, -0.34,  0.20,  0.26,  0.01, -0.07, -0.28, -0.33,  0.09,  0.00],
 [-0.06,  0.17, -0.14, -0.31, -0.02, -0.14, -0.26, -0.16, -0.05,  0.21,  0.25, -0.13,  0.14,  0.29,  0.01, -0.20, -0.08,  0.34,  0.09,  0.18,  0.00],
 [-0.09, -0.35, -0.11, -0.29,  0.19, -0.14,  0.00, -0.26, -0.19,  0.14,  0.20, -0.09,  0.19,  0.31, -0.07, -0.08,  0.03,  0.22,  0.13,  0.25,  0.00],
 [-0.09, -0.16,  0.06,  0.24,  0.08,  0.08,  0.29,  0.18, -0.12,  0.02, -0.09,  0.22, -0.67, -0.16, -0.28,  0.34,  0.22, -0.12, -0.04, -0.07,  0.00],
 [ 0.09, -0.25, -0.20,  0.00,  0.04, -0.20, -0.10,  0.14, -0.34,  0.11,  0.24, -0.21, -0.13,  0.00, -0.33,  0.09,  0.13, -0.04, -0.06,  0.02,  0.00],
 [-0.10,  0.30,  0.50,  0.58,  0.06,  0.24,  0.34,  0.16,  0.19, -0.25, -0.29,  0.44, -0.14, -0.22,  0.09,  0.18,  0.25, -0.07,  0.02, -0.29,  0.00],
 [ 0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00]])

#kt = (sc.R/4184)*(20+273.15)
#import scipy.constants as sc
kt = 0.6


############# Useful Functions #############
def aa_to_cdn_seq(aa_seq):
    '''
    Given an amino acid sequence, returns a codon sequence by randomly drawing from the syn codons
    '''
    cdn_seq = []
    for aa in aa_seq:
        aa_idx = AminoAcid[aa]
        syn_cdns = AA_SYNCDNS[AminoAcid[aa]]
        cdn_idx = syn_cdns[np.random.randint(0,len(syn_cdns))]
        cdn_seq.append(Codon[cdn_idx])
    return("".join(cdn_seq))
    
def aa_index_to_aa(aa_idx_seq):
    '''
        Given an aa idx seqeucne, convert it to a aa sequence. 
    '''
    aa_seq = []
    for i in range(len(aa_idx_seq)):
        aa = aa_idx_seq[i]
        for idx in AminoAcid:
            if AminoAcid[idx] == aa:
                aa_seq.append(idx)
    return "".join(aa_seq)

	
def Codon_to_AA(seqCodon):
    '''
        Given a codon index sequence, convert it to an amino acid index sequence.
    '''
    return [ AminoAcid[Codon_AA[Codon[cdn]]] for cdn in seqCodon ]

def Codon_to_Cindex(CodonSeq):
    '''
        Given a codon sequence, convert it to a codon index sequence. 
    '''
    seqCindex = []
    for i in range(0,len(CodonSeq),3):
        cdn = CodonSeq[i:i+3]
        for idx in Codon:
            if Codon[idx] == cdn:
                seqCindex.append(idx)
    return seqCindex

def AminoAcid_to_AAindex(Seq):
    '''
        Given an Amino acid sequence, convert it to an AA index sequence.
    '''
    SeqIndex = [ AminoAcid[i] for i in Seq]
    return (SeqIndex)

def NumberOfContacts(ContactMap):
    '''
        Calculates the number of contacts per site given a contact map (structure).
    '''
    Contacts = []
    for site in range(len(ContactMap)):
        Contacts.append(sum(ContactMap[site]))
    return Contacts        

def write_to_excel(file, LIST):
    '''
        writes a list to a csv file
    '''
    with open(file, 'w') as file:
        wr = csv.writer(file, quoting=csv.QUOTE_ALL)
        wr.writerow(LIST)
    print('Done')



############# Energy Calculations #############
@nb.jit(nopython = True)
def getG(aaSeq, ContactMap):
    '''
        Calculates the free energy of a sequence given a structure.
            Free energy is the sum of potentials between 
                amino acids in contact. 
        
        ContactMap: an array (len(seq)xlen(seq)) specifying residues in contact
        Sequence: Amino acid index sequence
    '''
    S = list(aaSeq)
    G = 0 
    for site1 in range(len(aaSeq)):
        for site2 in range(site1, len(aaSeq)):
            if ContactMap[site1,site2] == True:
                G +=  MJ85[S[site1], S[site2]]
    return G

def getGalt(aaSeq):
    '''
        Calculates the free energy of a given sequence in all alternative structures.
        
        Returns G_alt: list of free energies in each alt structure (ordered).
    '''
    Galt = [getG(aaSeq, AltMap) for AltMap in ContactMaps.values()]
    return Galt


@nb.jit
def Fitness(Sequence, isCodon, ContactMapNS):
    '''
        Sequence:index sequence
        isCodon: specify if index sequence is codon or AA; 1 => codon, 0 => AA
    '''
    SEQ = list(Sequence)
    if isCodon:
        SEQ = Codon_to_AA(SEQ)
    Gns = getG(SEQ, ContactMapNS) #G in Native structure
    Galt =  getGalt(SEQ) #G in alt structures
    dG = Gns - np.mean(Galt) + (np.var(Galt) / (2*kt)) +(kt*np.log(3.4**len(SEQ))) #delta G
    fitness = np.exp(-1.*dG/kt)/(1+ np.exp(-1.*dG/kt)) #Prob of folding = fitness
    return(fitness, Gns, Galt, dG)



@nb.jit(nopython = True)
def simple_getG(site, G, new_aaSeq, old_aaSeq, ContactMap):
    '''
        site: the site where newSeq and oldSeq differ
        G: free energy of oldSeq in ContactMap
        new_aaSeq, old_aaSeq: Amino acid idx sequence 
    '''
    Gnew = G
    for site2 in range(len(old_aaSeq)):
        if ContactMap[site,site2] == True:
            # subtract old interactions and add new ones 
            Gnew +=  MJ85[new_aaSeq[site],new_aaSeq[site2]]- MJ85[old_aaSeq[site],old_aaSeq[site2]] 
    return(Gnew)

def simple_getGalt(site, G, new_aaSeq, old_aaSeq):
    '''
        site: the site where newSeq and oldSeq differ
        G: free energy of oldSeq in ContactMap
        new_aaSeq, old_aaSeq: Amino acid idx sequence 
    '''
    Galt = [simple_getG(site, G[i], new_aaSeq, old_aaSeq, AltMap) for i, AltMap in enumerate(ContactMaps.values())]
    return Galt

@nb.jit
def simple_Fitness(newSeq, site, Gns, Galt, oldSeq, ContactMapNS):
    '''
        newSeq, oldSeq: codon idx sequences
        site: the site where newSeq and oldSeq differ
        G: free energy of oldSeq in ContactMap
        Galt: list of free energies in Alt Structures
    '''
    SEQ = list(newSeq); SEQ = Codon_to_AA(SEQ)
    SEQold = list(oldSeq); SEQold = Codon_to_AA(SEQold)
    Gns = simple_getG(site, Gns, SEQ, SEQold, ContactMapNS) #G in Native structure
    Galt =  simple_getGalt(site, Galt, SEQ, SEQold) #G in alt structures
    dG = Gns - np.mean(Galt) + (np.var(Galt) / (2*kt)) +(kt*np.log(3.4**len(SEQ))) #delta G
    fitness = np.exp(-1.*dG/kt)/(1+ np.exp(-1.*dG/kt)) #Prob of folding = fitness
    return(fitness)


############# EVOLUTIONARY MODEL #############
def MutationMatrix(pi_nuc, a,b,c,d,e,f):
    '''
        Returns mutation matrix. 
        
        pi_nuc = nucleotide frequencies <T,C,A,G>
        a,b,c,d,e,f = transition rates between nucelotides  
    '''
    GTR = np.dot(np.matrix([[0, a, b, c], 
                             [a, 0, d, e], 
                             [b, d, 0, f], 
                             [c, e, f, 0]]), np.diag(pi_nuc))
    #Fill in diagonal elements s.t. row sums to 0 
    for i in range(0, len(GTR)):
        GTR[i,i] = -np.sum(GTR[i,:])
    return GTR

def NucDiff(source, target):
    '''
        Get the nucleotide difference(s) between two codons.
        
        Returns a string of nucleotide difference between source and target codons
        source, target = three-letter sting represting a codons 
    '''
    return "".join( [source[i]+target[i] for i in range(len(source)) if source[i] != target[i]] )

def SingleStepNeighbour():
    '''
        makes a dictionary with codons as keys and thier single step neighbours as values
    '''
    SingleStepNeighbour = {}
    for codon in range(61):
        SingleStepNeighbour[codon] = []
        for codon2 in range(61):
            diff = NucDiff(Codon[codon], Codon[codon2])
            if len(diff) == 2:
                SingleStepNeighbour[codon].append(codon2)
    return (SingleStepNeighbour)
SingleStepNeighbour = SingleStepNeighbour()


def NonSyn_SingleStepNeighbour():
    '''
        makes a dictionary with codons as keys and 
        thier single step nonsynonymous neighbours as values
    '''
    NonSyn_SingleStepNeighbour = {}
    for codon in range(61):
        NonSyn_SingleStepNeighbour[codon] = []
        for codon2 in range(61):
            diff = NucDiff(Codon[codon], Codon[codon2])
            if len(diff) == 2 and Codon_AA[Codon[codon]] != Codon_AA[Codon[codon2]]:
                NonSyn_SingleStepNeighbour[codon].append(codon2)
    return (NonSyn_SingleStepNeighbour)

NonSynSingleStepNeighbour = NonSyn_SingleStepNeighbour()


def SynRateMatrix(GTR):
    '''
     Returns matrix with synonymous rates 
    '''
    SynM = np.zeros((61, 61))
    for codon in range(61):
        for codon2 in range(61):
            diff = NucDiff(Codon[codon], Codon[codon2])
            if len(diff) == 2 and Codon_AA[Codon[codon]] == Codon_AA[Codon[codon2]]:
                n1 = Nucleotide[diff[0]]; n2 = Nucleotide[diff[1]]
                SynM[codon, codon2] = GTR[n1,n2]
    return(SynM)

                
@nb.jit
def SiteSpecificTransitionMatrix(site, CurrentSeq, Gns, Galt, Neff, GTR, SynM, ContactMapNS):
    '''
        Creates site specific transition MATRIX Q^h (61x61)
       
        CurrentSeq: codon idx sequence
    '''
    Q = np.zeros((61,61))
    for CurrentCodon in range(61):
        for newCodon in NonSynSingleStepNeighbour[CurrentCodon]:
            diff = NucDiff(Codon[CurrentCodon], Codon[newCodon])
            n1 = Nucleotide[diff[0]]; n2 = Nucleotide[diff[1]]

            Seq = list(CurrentSeq); Seq[site] = CurrentCodon
            NewSeq = list(CurrentSeq); NewSeq[site] = newCodon
            Fit = simple_Fitness(Seq, site, Gns, list(Galt), CurrentSeq, ContactMapNS)
            newFit = simple_Fitness(NewSeq, site, Gns, list(Galt), CurrentSeq, ContactMapNS)
            Sij = 2*Neff*(newFit - Fit)
            #Sij = (newFit - Fit)/Fit #G&P way of doing it 
            if abs(Sij) <= 1e-10:
                Q[CurrentCodon, newCodon] = GTR[n1,n2]
            else:
                Q[CurrentCodon, newCodon] = GTR[n1,n2] * (Sij)/(1 - np.exp(-1.*Sij))
                #Q[CurrentCodon, newCodon] = GTR[n1,n2]* 2 * Neff* ( (1 - np.exp(-2*Sij))/(1-np.exp(-4*Neff*Sij))) #G&P way of doing it 
                
    return (Q + SynM)

def SiteSpecificTransitionMatrix_AllSites(CurrentSeq, num_cores, Neff, GTR, SynM, ContactMapNS):
    '''
        Creates a SiteSpecific transition matrix Q^h for all sites.
                
        CurrentSeq: Codon index Sequence
    '''

    Sold = list(CurrentSeq)
    Pold, GNSold, GALTold, dGold = Fitness(Sold, 1, ContactMapNS)
    results = []
    results = Parallel(n_jobs=num_cores)(delayed(SiteSpecificTransitionMatrix)(i, Sold, GNSold, GALTold, Neff, GTR, SynM, ContactMapNS) for i in range(0, len(Sold)) )
    Q = np.zeros((len(CurrentSeq),61,61))

    for Idx, i in enumerate(results):
        Q[Idx] = i
    return(Q)

      
@nb.jit    
def SiteSpecificTransitionVector_SS(site, CurrentSeq, Fit, Gns, Galt, Neff, GTR, ContactMapNS):
    '''
        Site specific transition vector
        Keeps track of Sijs 
       
        CurrentSeq: codon idx sequence
    '''
    Q = np.zeros(61)
    M = np.zeros(61)
    SS = []
    CurrentCodon = CurrentSeq[site]
    for newCodon in SingleStepNeighbour[CurrentCodon]:
        diff = NucDiff(Codon[CurrentCodon], Codon[newCodon])
        n1 = Nucleotide[diff[0]]; n2 = Nucleotide[diff[1]]
        M[newCodon] = GTR[n1,n2]
        if Codon_AA[Codon[newCodon]] == Codon_AA[Codon[CurrentCodon]]:
            Q[newCodon] = GTR[n1,n2]
            SS.append(0)
        else:
            NewSeq = list(CurrentSeq); NewSeq[site] = newCodon
            newFit = simple_Fitness(NewSeq, site, Gns, list(Galt), CurrentSeq, ContactMapNS)
#            Sij = (newFit - Fit)/Fit
            Sij = 2*Neff*(newFit - Fit)
            SS.append(Sij)
            if abs(Sij) <= 1e-10:
                Q[newCodon] = GTR[n1,n2]
            else:
                Q[newCodon] = GTR[n1,n2] * (Sij)/(1 - np.exp(-1.*Sij))
                #Q[newCodon] = GTR[n1,n2]* 2 * Neff* ( (1 - np.exp(-2*Sij))/(1-np.exp(-4*Neff*Sij)))
    return ((Q,M, SS))


def TransitionVector_SS(CurrentSeq, num_cores, Neff, GTR, ContactMapNS):
    '''
        Creates a row of transition matrix.
        Keeps track of Sijs 
                
        CurrentSeq: Codon index Sequence
    '''

    Sold = list(CurrentSeq)
    Pold, GNSold, GALTold, dGold = Fitness(Sold, 1, ContactMapNS)
    results = []
    results = Parallel(n_jobs=num_cores)(delayed(SiteSpecificTransitionVector_SS)(i, Sold,  Pold, GNSold, GALTold, Neff, GTR, ContactMapNS) for i in range(0, len(Sold)) )
    Q = np.zeros((len(CurrentSeq),61))
    M = np.zeros((len(CurrentSeq),61))
    SIJ = []

    for Idx, i in enumerate(results):
        Q[Idx] = i[0]
        M[Idx] = i[1]
        SIJ.append(i[2])
    return(Q, M, SIJ, Pold, dGold)
    
@nb.jit    
def SiteSpecificTransitionVector(site, CurrentSeq, Fit, Gns, Galt, Neff, GTR, ContactMapNS):
    '''
        Site specific transition vector
        
        CurrentSeq: codon idx sequence
    '''
    Q = np.zeros(61)
    M = np.zeros(61)
    CurrentCodon = CurrentSeq[site]
    for newCodon in SingleStepNeighbour[CurrentCodon]:
        diff = NucDiff(Codon[CurrentCodon], Codon[newCodon])
        n1 = Nucleotide[diff[0]]; n2 = Nucleotide[diff[1]]
        M[newCodon] = GTR[n1,n2]
        if Codon_AA[Codon[newCodon]] == Codon_AA[Codon[CurrentCodon]]:
            Q[newCodon] = GTR[n1,n2]
        else:
            NewSeq = list(CurrentSeq); NewSeq[site] = newCodon
            newFit = simple_Fitness(NewSeq, site, Gns, list(Galt), CurrentSeq, ContactMapNS)
#            Sij = (newFit - Fit)/Fit
            Sij = 2*Neff*(newFit - Fit)
            if abs(Sij) <= 1e-10:
                Q[newCodon] = GTR[n1,n2]
            else:
                Q[newCodon] = GTR[n1,n2] * (Sij)/(1 - np.exp(-1.*Sij))
#                Q[newCodon] = GTR[n1,n2]* 2 * Neff* ( (1 - np.exp(-2*Sij))/(1-np.exp(-4*Neff*Sij)))
    return ((Q,M))


def TransitionVector(CurrentSeq, num_cores, Neff, GTR, ContactMapNS):
    '''
        Creates a row of transition matrix.
                
        CurrentSeq: Codon index Sequence
    '''

    Sold = list(CurrentSeq)
    Pold, GNSold, GALTold, dGold = Fitness(Sold, 1, ContactMapNS)
    results = []
    results = Parallel(n_jobs=num_cores)(delayed(SiteSpecificTransitionVector)(i, Sold,  Pold, GNSold, GALTold, Neff, GTR, ContactMapNS) for i in range(0, len(Sold)) )
    Q = np.zeros((len(CurrentSeq),61))
    M = np.zeros((len(CurrentSeq),61))

    for Idx, i in enumerate(results):
        Q[Idx] = i[0]
        M[Idx] = i[1]
    return(Q)


def IndicatorMatrix(Sold):
    IN = np.zeros((len(Sold), 61))
    IS = np.zeros((len(Sold), 61))
    for site in range(len(Sold)):
        for newCodon in range(61):
            diff = NucDiff(Codon[Sold[site]], Codon[newCodon])
            if len(diff) != 2:
                pass
            else:
                oldAA = Codon_AA[Codon[Sold[site]]]
                newAA = Codon_AA[Codon[newCodon]]
                if oldAA == newAA:
                    IS[site, newCodon] = 1
                else:
                    IN[site, newCodon] = 1
    return IN, IS
                
    
@nb.jit(nb.types.Tuple((nb.int_, nb.int_, nb.float64))(nb.float64[:,:]))
def NextSubstitution(Q):
    rate = np.sum(Q)
    TransitionProb = Q / rate
    X = np.random.multinomial(1, TransitionProb.reshape(len(TransitionProb)* 61)).reshape(len(TransitionProb), 61)
    site, codon = np.unravel_index(X.argmax(), X.shape)
    return (site, codon, rate)

@nb.jit
def NextSubstitution_FullQh(Q, Seq):
    '''
        Given the full Q (len(Seq)X61x61) matrix calculates the most likely site and codon
    '''
    #Get Qv the vectorized transition matrix (row of Q for currentcodon)
    Qv = np.zeros((len(Seq),61))
    for site in range(len(Seq)):
        CurrentCodon = Seq[site]
        Qv[site] = Q[site, CurrentCodon]
    #Site and Codon change
    rate = np.sum(Qv)
    TransitionProb = Qv / rate
    X = np.random.multinomial(1, TransitionProb.reshape(len(TransitionProb)* 61)).reshape(len(TransitionProb), 61)
    site, codon = np.unravel_index(X.argmax(), X.shape)
    return (site, codon)


#SequenceFitness = {}
#def CheckFit(Sequence):
#    '''
#        Checks if fitness calculation was already preformed for given Sequence
#    '''
#    if tuple(Sequence) in SequenceFitness:
#        Pnew = SequenceFitness[tuple(Sequence)]
#    else:
#        Pnew = Fitness(Sequence, 1)[0]
#        SequenceFitness[tuple(Sequence)] = Pnew
#    return(Pnew)  
