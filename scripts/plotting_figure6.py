#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 13:30:00 2020

@author: Noory333
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 10:37:38 2019

@author: Noor Youssef

Plots the relative error in dN/dS versus omega for M3 k = 2
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Arial"
plt.rcParams["mathtext.fontset"] = 'cm'

    
#expeced site-specific  rates# 
def calculate_dnds_h(generating_model, protein):
    if generating_model == "Ne2_C-SI" or generating_model == "Ne2_S-SI":        
        Kn = np.genfromtxt(RstPath + 'dNdS/' +  protein + "_" + generating_model +'_Kn.csv', delimiter=' ')
        Ln = np.genfromtxt(RstPath + 'dNdS/' + protein + "_" + generating_model +'_Ln.csv', delimiter=' ')
        return(Kn, Ln)
    elif generating_model == "Ne2_S-SD":
        Kn = []; Ln = []
        for trial in range(1,51):
            full_Kn = np.genfromtxt(RstPath +'dNdS/' +  protein + "_" + generating_model + '_Kn_Ln/Kn_' + str(trial) + '.csv', delimiter=' ')
            full_Ln = np.genfromtxt(RstPath +'dNdS/' +  protein + "_" + generating_model + '_Kn_Ln/Ln_' + str(trial) + '.csv', delimiter=' ')
            Kn.append(np.sum(full_Kn, axis = 0))
            Ln.append(np.sum(full_Ln, axis = 0))
        return(np.array(Kn), np.array(Ln))


def expected_dNdS_M3k2(post_prob_1, post_prob_2, Kn, Ln):
    dnds_w1 = []; dnds_w2 = []
    for trial in range(50):
        kn_1 = [pp1*kn for pp1, kn in zip(post_prob_1[trial], Kn[trial])]
        ln_1 = [pp1*ln for pp1, ln in zip(post_prob_1[trial], Ln[trial])]
        dnds_w1.append(np.sum(kn_1) / np.sum(ln_1))
        
        kn_2 = [pp2*kn for pp2, kn in zip(post_prob_2[trial], Kn[trial])]
        ln_2 = [pp2*ln for pp2, ln in zip(post_prob_2[trial], Ln[trial])]
        dnds_w2.append(np.sum(kn_2) / np.sum(ln_2))
    return(np.array(dnds_w1), np.array(dnds_w2))


def _main_results_(protein, gen_model):
    #omega# 
    M3k2 = pd.read_csv(RstPath + "Inference/"+ protein + "_" + gen_model +"_M3_k2.csv", sep = " ")
    
    #dN/dS_1 dN/dS_2 # 
    full_Kn, full_Ln = calculate_dnds_h(gen_model, protein)
    
    post_prob_1 = np.genfromtxt(RstPath + 'post_prob_1_'+ protein + "_" + gen_model +'_' + 'M3_k2' + '.csv', delimiter=' ')
    post_prob_2 = np.genfromtxt(RstPath + 'post_prob_2_'+ protein + "_" + gen_model +'_' + 'M3_k2' + '.csv', delimiter=' ')
    dnds_1, dnds_2 = expected_dNdS_M3k2(post_prob_1, post_prob_2, full_Kn, full_Ln)

    # relative error
    rel_error_w1 = (M3k2.w0 / dnds_1) - 1
    rel_error_w2 = (M3k2.w1 / dnds_2) - 1

    return(rel_error_w1, rel_error_w2)



RstPath = '../Results/'

qhw_CSI_rel_error_w1, qhw_CSI_rel_error_w2 = _main_results_("1qhw", "Ne2_C-SI")
qhw_SSI_rel_error_w1, qhw_SSI_rel_error_w2 = _main_results_("1qhw", "Ne2_S-SI")
qhw_SSD_rel_error_w1, qhw_SSD_rel_error_w2 = _main_results_("1qhw", "Ne2_S-SD")

ppn_CSI_rel_error_w1, ppn_CSI_rel_error_w2 = _main_results_("2ppn", "Ne2_C-SI")
ppn_SSI_rel_error_w1, ppn_SSI_rel_error_w2 = _main_results_("2ppn", "Ne2_S-SI")
ppn_SSD_rel_error_w1, ppn_SSD_rel_error_w2 = _main_results_("2ppn", "Ne2_S-SD")

pek_CSI_rel_error_w1, pek_CSI_rel_error_w2 = _main_results_("1pek", "Ne2_C-SI")
pek_SSI_rel_error_w1, pek_SSI_rel_error_w2 = _main_results_("1pek", "Ne2_S-SI")
pek_SSD_rel_error_w1, pek_SSD_rel_error_w2 = _main_results_("1pek", "Ne2_S-SD")


#%%###########################################################
############               Plot                 ############
############################################################
from matplotlib.lines import Line2D
f, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey = True, sharex = True, figsize = (12, 4))

c = ["darkblue", "darkred", "chocolate"]

#1QHW# 
#omega_1# 
for d, (i,model) in enumerate(zip([1,2.5,4], [qhw_CSI_rel_error_w1, qhw_SSI_rel_error_w1, qhw_SSD_rel_error_w1])):
    y = list(model)
    x = np.random.normal(i, 0.04, size=len(y))
    ax1.errorbar(i,np.mean(y), np.std(y), color = c[d], marker = 'o', label = r"$\omega_1$")
    ax1.plot(x, y, 'o', markersize = 4, color = c[d], alpha=0.1)
 
#omega_2# 
for d, (i,model) in enumerate(zip([1.5,3,4.5], [qhw_CSI_rel_error_w2, qhw_SSI_rel_error_w2, qhw_SSD_rel_error_w2])):
    y = model
    x = np.random.normal(i, 0.04, size=len(y))
    ax1.errorbar(i,np.mean(y), np.std(y), color = c[d], marker = '^', label = r"$\omega_2$")
    ax1.plot(x, y, '^', markersize = 4, color = c[d], alpha=0.1)
ax1.axhline(0, linestyle = "--", c = 'k')
ax1.set_xticks([1.25,2.75,4.25])
ax1.set_xticklabels(["C-SI", "S-SI", "S-SD"], fontsize = 12)
ax1.set_title("1QHW", fontsize = 14)
ax1.set_ylabel("Relative error", fontsize = 12)

#2ppn# 
#omega_1# 
for d, (i,model) in enumerate(zip([1,2.5,4], [ppn_CSI_rel_error_w1, ppn_SSI_rel_error_w1, ppn_SSD_rel_error_w1])):
    y = list(model)
    x = np.random.normal(i, 0.04, size=len(y))
    ax2.errorbar(i,np.mean(y), np.std(y), color = c[d], marker = 'o', label = r"$\omega_1$")
    ax2.plot(x, y, 'o', markersize = 4, color = c[d], alpha=0.1)
 
#omega_2# 
for d, (i,model) in enumerate(zip([1.5,3,4.5], [ppn_CSI_rel_error_w2, ppn_SSI_rel_error_w2, ppn_SSD_rel_error_w2])):
    y = model
    x = np.random.normal(i, 0.04, size=len(y))
    ax2.errorbar(i,np.mean(y), np.std(y), color = c[d], marker = '^', label = r"$\omega_2$")
    ax2.plot(x, y, '^', markersize = 4, color = c[d], alpha=0.1)
ax2.axhline(0, linestyle = "--", c = 'k')
ax2.set_xticklabels(["C-SI", "S-SI", "S-SD"], fontsize = 12)
ax2.set_title("2PPN", fontsize = 14)

#1PEK# 
#omega_1# 
for d, (i,model) in enumerate(zip([1,2.5,4], [pek_CSI_rel_error_w1, pek_SSI_rel_error_w1, pek_SSD_rel_error_w1])):
    y = list(model)
    x = np.random.normal(i, 0.04, size=len(y))
    ax3.errorbar(i,np.mean(y), np.std(y), color = c[d], marker = 'o', label = r"$\omega_1$")
    ax3.plot(x, y, 'o', markersize = 4, color = c[d], alpha=0.1)
 
#omega_2# 
for d, (i,model) in enumerate(zip([1.5,3,4.5], [pek_CSI_rel_error_w2, pek_SSI_rel_error_w2, pek_SSD_rel_error_w2])):
    y = model
    x = np.random.normal(i, 0.04, size=len(y))
    ax3.errorbar(i,np.mean(y), np.std(y), color = c[d], marker = '^', label = r"$\omega_2$")
    ax3.plot(x, y, '^', markersize = 4, color = c[d], alpha=0.1)
ax3.axhline(0, linestyle = "--", c = 'k')
ax3.set_title("1PEK", fontsize = 14)
ax3.set_xticklabels(["C-SI", "S-SI", "S-SD"], fontsize = 12)
legend_elements = [Line2D([0], [0], marker='o', color='k', label=r"$\omega_1$"), Line2D([0], [0], marker='^', color='k', label=r"$\omega_2$")]
ax3.legend(handles=legend_elements)


plt.savefig('../Figures/Figure6.png', dpi = 450, bbox_inches = 'tight')



