#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 13:04:05 2020

@author: nooryoussef

Correlation between dNdS and wh 
"""
import numpy as np 
from scipy import stats

RstPath = '../Results/'

def calculate_dnds_h(generating_model, protein):
    '''
        Extracts expeced site-specific rates
    '''
    if generating_model == "Ne2_C-SI" or generating_model == "Ne2_S-SI":        
        Kn = np.genfromtxt(RstPath + 'dNdS/' +  protein + "_" + 
                           generating_model +'_Kn.csv', delimiter=' ')
        Ln = np.genfromtxt(RstPath + 'dNdS/' + protein + "_" + 
                           generating_model +'_Ln.csv', delimiter=' ')
        return(Kn, Ln)
    elif generating_model == "Ne2_S-SD":
        Kn = []; Ln = []
        for trial in range(1,51):
            full_Kn = np.genfromtxt(RstPath +'dNdS/' +  protein + "_" + 
                                    generating_model + '_Kn_Ln/Kn_' + 
                                    str(trial) + '.csv', delimiter=' ')
            full_Ln = np.genfromtxt(RstPath +'dNdS/' +  protein + "_" + 
                                    generating_model + '_Kn_Ln/Ln_' + 
                                    str(trial) + '.csv', delimiter=' ')
            Kn.append(np.sum(full_Kn, axis = 0))
            Ln.append(np.sum(full_Ln, axis = 0))
        return(np.array(Kn), np.array(Ln))
    
def correlation(dnds, omega):
    R = []; P = []
    for trial in range(50):
        r, p = stats.pearsonr(omega[trial], dnds[trial])
        R.append(r)
        P.append(p)
    return(np.array(R), np.array(P))
    

dnds = {}
omega = {}
for protein in ['1qhw', '2ppn', '1pek']:
    print('\n', protein, '\n')
    for model in ['Ne2_C-SI', 'Ne2_S-SI', 'Ne2_S-SD']:

        if protein == '2ppn':
            omega_h =  np.genfromtxt(RstPath + 'Inference/post_mean_omega_'+ protein + "_" + model + '_M3_k2.csv')
        else:
            omega_h =  np.genfromtxt(RstPath + 'Inference/post_mean_omega_'+ protein + "_" + model + '_M3_k3.csv')
        
        omega[protein + "_" + model] = omega_h
        
        #dN^h/dS^h#
        full_Kn, full_Ln = calculate_dnds_h(model, protein)
        dnds_h = full_Kn / full_Ln
        dnds[protein + "_" + model] = dnds_h
        
        R, P = correlation(dnds_h, omega_h)
        sig = [P <= 0.005]

        print (str(model) + " \n mean R (std) Pvalue " + 
        "%.3f" % np.mean(R) + "(" + "%.3f" % np.std(R) + ") " + 
        "%.3f" % np.mean(P) + " \n Number of significant trials " +
        str(np.sum(sig)))


        
