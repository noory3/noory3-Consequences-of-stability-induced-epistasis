#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 09:34:06 2020

@author: nooryoussef

Plotting distributions of dN^h/dS^h and omega MLEs
"""
import numpy as np
import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Arial"
plt.rcParams["mathtext.fontset"] = 'cm'

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
    
def _main_results_(protein, gen_model, M_model):
    # dNdS^h #
    full_Kn, full_Ln = calculate_dnds_h(gen_model, protein)
    dnds_h = full_Kn / full_Ln
    
    # Best fitting M-series model
    Mseries = pd.read_csv(RstPath + "Inference/" + protein + "_" + gen_model +"_"+ M_model+ ".csv", sep = " ")
    return (dnds_h, Mseries)

#%%###############################
#1qhw#
dnds_h_1qhw_CSI, Mmodel_1qhw_CSI = _main_results_('1qhw', "Ne2_C-SI", "M3_k2")
dnds_h_1qhw_SSI, Mmodel_1qhw_SSI = _main_results_('1qhw', "Ne2_S-SI", "M3_k2")
dnds_h_1qhw_SSD, Mmodel_1qhw_SSD = _main_results_('1qhw', "Ne2_S-SD", "M3_k2")
 
#2ppn# 
dnds_h_2ppn_CSI, Mmodel_2ppn_CSI = _main_results_('2ppn', "Ne2_C-SI", "M3_k2")
dnds_h_2ppn_SSI, Mmodel_2ppn_SSI = _main_results_('2ppn', "Ne2_S-SI", "M3_k2")
dnds_h_2ppn_SSD, Mmodel_2ppn_SSD = _main_results_('2ppn', "Ne2_S-SD", "M3_k2")
 
#1pek# 
dnds_h_1pek_CSI, Mmodel_1pek_CSI = _main_results_('1pek', "Ne2_C-SI", "M3_k3")
dnds_h_1pek_SSI, Mmodel_1pek_SSI = _main_results_('1pek', "Ne2_S-SI", "M3_k3")
dnds_h_1pek_SSD, Mmodel_1pek_SSD = _main_results_('1pek', "Ne2_S-SD", "M3_k3")

#%%######## plot ######### 
# C-SI Plots 
qhw = [dnds_h_1qhw_CSI, Mmodel_1qhw_CSI]
ppn = [dnds_h_2ppn_CSI, Mmodel_2ppn_CSI]
pek = [dnds_h_1pek_CSI, Mmodel_1pek_CSI]

f, ((ax4, ax5, ax6), (ax1, ax2 , ax3)) = plt.subplots(2, 3, sharex = True, sharey = False, gridspec_kw={'height_ratios': [0.7, 1]}, figsize = (8,2))

f.subplots_adjust(wspace=0.2, hspace = -0.2)

# 1qhw # 
sns.distplot(qhw[0].flatten(), ax = ax1, color = 'darkblue', kde=False, norm_hist=True, bins = 30)
ax1.set_yticklabels([])
ax1.set_ylabel('C-SI        ', fontsize = 16, rotation = 360, color = 'darkblue')
ax1.set_frame_on(False)
ax1.axes.get_xaxis().set_visible(True)
ax1.axhline(0, color = 'k', linewidth = 3)

ax1.set_yticks([])

box1 = ax4.boxplot([qhw[1].w0, qhw[1].w1], 
             positions=[1,1], showfliers=False, vert=False)
for item in ['boxes', 'whiskers', 'fliers', 'medians', 'caps']:
        plt.setp(box1[item], color='darkblue')

ax4.set_frame_on(False)
ax4.set_yticks([])
ax4.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
ax4.set_title('1QHW', fontsize = 16)

# 2ppn # 
sns.distplot(ppn[0].flatten(), ax = ax2, color = 'darkblue', kde=False, norm_hist=True, bins = 30)
ax2.set_yticks([])
ax2.set_frame_on(False)
ax2.axhline(0, color = 'k', linewidth = 3)

box1 = ax5.boxplot([ppn[1].w0, ppn[1].w1], 
             positions=[1,1], showfliers=False, vert=False)
for item in ['boxes', 'whiskers', 'fliers', 'medians', 'caps']:
        plt.setp(box1[item], color='darkblue')
ax5.set_yticks([])
ax5.set_frame_on(False)
ax5.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
ax5.set_title('2PPN', fontsize = 16)

# 1Pek # 
sns.distplot(pek[0].flatten(), ax = ax3, color = 'darkblue', kde=False, norm_hist=True,bins = 30)
ax3.set_yticks([])
ax3.legend(fontsize = 8)
ax3.set_frame_on(False)
ax3.axhline(0, color = 'k', linewidth = 3)

box1 = ax6.boxplot([pek[1].w0, pek[1].w1, pek[1].w2], 
             positions=[1,1,1], showfliers=False, vert=False)
for item in ['boxes', 'whiskers', 'fliers', 'medians', 'caps']:
        plt.setp(box1[item], color='darkblue')
ax6.set_yticks([])
ax6.set_frame_on(False)
ax6.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
ax6.set_title('1PEK', fontsize = 16)

ax6.set_xticks([0,0.5,1, 1.5])
ax6.set_xlim(-0.05,1.7)

plt.legend(frameon=False)
plt.savefig('../Figures/C-SI_dNdSh_mle.png', dpi = 450, bbox_inches = 'tight')

#%%# S-SI Plots 
qhw = [dnds_h_1qhw_SSI, Mmodel_1qhw_SSI]
ppn = [dnds_h_2ppn_SSI, Mmodel_2ppn_SSI]
pek = [dnds_h_1pek_SSI, Mmodel_1pek_SSI]

f, ((ax4, ax5, ax6), (ax1, ax2 , ax3)) = plt.subplots(2, 3, sharex = True, sharey = False, gridspec_kw={'height_ratios': [0.7, 1]}, figsize = (8,2))

f.subplots_adjust(wspace=0.2, hspace = -0.2)

# 1qhw # 
sns.distplot(qhw[0].flatten(), ax = ax1, color = 'darkred', kde=False, norm_hist=True, bins = 30)
ax1.set_yticklabels([])
ax1.set_ylabel('S-SI        ', fontsize = 16, rotation = 360, color = 'darkred')
ax1.set_frame_on(False)
ax1.axes.get_xaxis().set_visible(True)
ax1.axhline(0, color = 'k', linewidth = 3)

ax1.set_yticks([])
box1 = ax4.boxplot([qhw[1].w0, qhw[1].w1], 
             positions=[1,1], showfliers=False, vert=False)
for item in ['boxes', 'whiskers', 'fliers', 'medians', 'caps']:
        plt.setp(box1[item], color='darkred')

ax4.set_frame_on(False)
ax4.set_yticks([])
ax4.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)


# 2ppn # 
sns.distplot(ppn[0].flatten(), ax = ax2, color = 'darkred', kde=False, norm_hist=True, bins = 30)
ax2.set_yticks([])
ax2.set_frame_on(False)
ax2.axhline(0, color = 'k', linewidth = 3)

box1 = ax5.boxplot([ppn[1].w0, ppn[1].w1], 
             positions=[1,1], showfliers=False, vert=False)
for item in ['boxes', 'whiskers', 'fliers', 'medians', 'caps']:
        plt.setp(box1[item], color='darkred')
ax5.set_yticks([])
ax5.set_frame_on(False)
ax5.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)

# 1Pek # 
sns.distplot(pek[0].flatten(), ax = ax3, color = 'darkred', kde=False, norm_hist=True,bins = 30)
ax3.set_yticks([])
ax3.legend(fontsize = 8)
ax3.set_frame_on(False)
ax3.axhline(0, color = 'k', linewidth = 3)

box1 = ax6.boxplot([pek[1].w0, pek[1].w1, pek[1].w2], 
             positions=[1,1,1], showfliers=False, vert=False)
for item in ['boxes', 'whiskers', 'fliers', 'medians', 'caps']:
        plt.setp(box1[item], color='darkred')
ax6.set_yticks([])
ax6.set_frame_on(False)
ax6.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)

ax6.set_xticks([0,0.5,1])
ax6.set_xlim(-0.05,1)
plt.legend(frameon=False)

plt.savefig('../Figures/S-SI_dndsh_mles.png', dpi = 450, bbox_inches = 'tight')

#%% S-SD Plots
qhw = [dnds_h_1qhw_SSD, Mmodel_1qhw_SSD]
ppn = [dnds_h_2ppn_SSD, Mmodel_2ppn_SSD]
pek = [dnds_h_1pek_SSD, Mmodel_1pek_SSD]

f, ((ax4, ax5, ax6), (ax1, ax2 , ax3)) = plt.subplots(2, 3, sharex = True, sharey = False, gridspec_kw={'height_ratios': [0.7, 1]}, figsize = (8,2))

f.subplots_adjust(wspace=0.2, hspace = -0.2)

# 1qhw # 
sns.distplot(qhw[0].flatten(), ax = ax1, color = 'chocolate', kde=False, norm_hist=True, bins = 30)
ax1.set_yticklabels([])
ax1.set_ylabel('S-SD        ', fontsize = 16, rotation = 360, color = 'chocolate')
ax1.set_xlabel('$dN^h/dS^h$')
ax1.set_frame_on(False)
ax1.axes.get_xaxis().set_visible(True)
ax1.axhline(0, color = 'k', linewidth = 3)

ax1.set_yticks([])

box1 = ax4.boxplot([qhw[1].w0, qhw[1].w1], 
             positions=[1,1], showfliers=False, vert=False)
for item in ['boxes', 'whiskers', 'fliers', 'medians', 'caps']:
        plt.setp(box1[item], color='chocolate')
ax4.set_frame_on(False)
ax4.set_yticks([])
ax4.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)


# 2ppn # 
sns.distplot(ppn[0].flatten(), ax = ax2, color = 'chocolate', kde=False, norm_hist=True, bins = 30)
ax2.set_yticks([])
ax2.set_xlabel('$dN^h/dS^h$')
ax2.set_frame_on(False)
ax2.axhline(0, color = 'k', linewidth = 3)

box1 = ax5.boxplot([ppn[1].w0, ppn[1].w1], 
             positions=[1,1], showfliers=False, vert=False)
for item in ['boxes', 'whiskers', 'fliers', 'medians', 'caps']:
        plt.setp(box1[item], color='chocolate')
ax5.set_yticks([])
ax5.set_frame_on(False)
ax5.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)

# 1Pek # 
sns.distplot(pek[0].flatten(), ax = ax3, color = 'chocolate', kde=False, norm_hist=True,bins = 30)
ax3.set_yticks([])
ax3.set_xlabel('$dN^h/dS^h$')
ax3.legend(fontsize = 8)
ax3.set_frame_on(False)
ax3.axhline(0, color = 'k', linewidth = 3)

box1 = ax6.boxplot([pek[1].w0, pek[1].w1, pek[1].w2], 
             positions=[1,1,1], showfliers=False, vert=False)
for item in ['boxes', 'whiskers', 'fliers', 'medians', 'caps']:
        plt.setp(box1[item], color='chocolate')
ax6.set_yticks([])
ax6.set_frame_on(False)
ax6.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)

ax6.set_xticks([0,0.5,1])
ax6.set_xlim(-0.05,1)
plt.legend(frameon=False)

plt.savefig('../Figures/S-SD_dndsh_mle.png', dpi = 450, bbox_inches = 'tight')
