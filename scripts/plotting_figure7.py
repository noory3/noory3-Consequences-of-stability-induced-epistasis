#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 10:41:58 2020

@author: nooryoussef

Plotting dN^h/dS^h S-SD versus dN^h/dS^h S-SI 
And rate difference 
"""
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Arial"
plt.rcParams["mathtext.fontset"] = 'cm'

RstPath = '../Results/dNdS/'

#### Expected Rates ####
def calculate_dnds_h(generating_model, protein):
    if generating_model == "Ne2_C-SI" or generating_model == "Ne2_S-SI":        
        Kn = np.genfromtxt(RstPath + protein + "_" +  generating_model +'_Kn.csv', delimiter=' ')
        Ln = np.genfromtxt(RstPath + protein + "_" +  generating_model +'_Ln.csv', delimiter=' ')
        return(Kn, Ln)
    elif generating_model == "Ne2_S-SD":
        Kn = []; Ln = []
        for trial in range(1,51):
            full_Kn = np.genfromtxt(RstPath + protein + "_" + generating_model + '_Kn_Ln/Kn_' + str(trial) + '.csv', delimiter=' ')
            full_Ln = np.genfromtxt(RstPath + protein + "_" + generating_model + '_Kn_Ln/Ln_' + str(trial) + '.csv', delimiter=' ')
            Kn.append(np.sum(full_Kn, axis = 0))
            Ln.append(np.sum(full_Ln, axis = 0))
        return(np.array(Kn), np.array(Ln))

def _main_results_(protein, gen_model):
    # dNdS^h #
    full_Kn, full_Ln = calculate_dnds_h(gen_model, protein)
    dnds_h = full_Kn / full_Ln
    return (dnds_h)

def create_bin(lower_bound, width, quantity):
    bins = []
    for low in np.arange(lower_bound, 
                     lower_bound + quantity*width, width):
        bins.append((low, low+width))
    return bins

def find_bin(value, bins):
    """ bins is a list of tuples, like [(0,20), (20, 40), (40, 60)],
        binning returns the smallest index i of bins so that
        bin[i][0] <= value < bin[i][1]
    """
    
    for i in range(0, len(bins)):
        if bins[i][0] <= value < bins[i][1]:
            return i
    return -1


def get_mean_and_std(dnds_SSI, dnds_SSD):
    d = np.zeros((50, dnds_SSI.shape[1]))
    for trial in range(50):
        dd = [dnds_SSD[trial][x] - dnds_SSI[trial][x] for x in range(dnds_SSI.shape[1])]
        d[trial] = dd
        
    bins = create_bin(0, 0.025, 41)
    binned_delta_dnds = [[] for x in range(len(bins))]
    binned_SSI_dnds = [[] for x in range(len(bins))]
    
    for delta_dnds, SSI_dnds in zip(d.flatten(),dnds_SSI.flatten()):
        bin_index = find_bin(SSI_dnds, bins)
        binned_SSI_dnds[bin_index].append(SSI_dnds)
        binned_delta_dnds[bin_index].append(delta_dnds)
    return(binned_SSI_dnds, binned_delta_dnds)

#%%###############################
#1qhw#
dnds_h_1qhw_SSI = _main_results_('1qhw', "Ne2_S-SI")
dnds_h_1qhw_SSD = _main_results_('1qhw', "Ne2_S-SD")
qhw_binned_SSI_dnds, qhw_binned_delta_dnds = get_mean_and_std(dnds_h_1qhw_SSI, dnds_h_1qhw_SSD)

#2ppn# 
dnds_h_2ppn_SSI = _main_results_('2ppn', "Ne2_S-SI")
dnds_h_2ppn_SSD = _main_results_('2ppn', "Ne2_S-SD")
ppn_binned_SSI_dnds, ppn_binned_delta_dnds = get_mean_and_std(dnds_h_2ppn_SSI, dnds_h_2ppn_SSD)


#1pek# 
dnds_h_1pek_SSI = _main_results_('1pek', "Ne2_S-SI")
dnds_h_1pek_SSD = _main_results_('1pek', "Ne2_S-SD")
pek_binned_SSI_dnds, pek_binned_delta_dnds = get_mean_and_std(dnds_h_1pek_SSI, dnds_h_1pek_SSD)


#%%
line_x = np.linspace(- 0.05, 1.05, 1000)

fig, ((ax1, ax2 , ax3), (ax4, ax5, ax6)) = plt.subplots(2,3, figsize = (15,6), sharex = 'col', sharey = False)

## 1QHW  ##
for trial in range(50):
    ax1.scatter(dnds_h_1qhw_SSI[trial], dnds_h_1qhw_SSD[trial],color = '0.5', alpha = 0.1, s = 10)

ax1.scatter(dnds_h_1qhw_SSI[0][2], dnds_h_1qhw_SSD[0][2],color = '0.5', alpha = 0.5, s = 30, label = r'$dN^h/dS^h$')
ax1.set_ylabel(r'expected $dN^h/dS^h$ under S-SD' , fontsize = 12)
ax1.plot(line_x, line_x, linestyle = '--', color = 'k', marker = '' , label = 'x = y')
ax1.set_ylim(-0.05, 1.05)
ax1.set_xlim(-0.05, 1.05)
ax1.set_title('1QHW', fontsize = 16)


ax4.set_ylabel('rate difference', fontsize = 12)  
for x, y in zip(qhw_binned_SSI_dnds, qhw_binned_delta_dnds):
    ax4.errorbar(np.mean(x),np.mean(y), yerr = np.std(y), xerr = np.std(x), color = 'k', marker = 'o')
ax4.axhline(linestyle = '--', color = 'k')
ax4.set_xlabel(r'expected $dN^h/dS^h$ under S-SI' , fontsize = 12)
ax4.set_yticks(np.arange(-0.2, 0.4, 0.1))

## 2PPN  ##
for trial in range(50):
    ax2.scatter(dnds_h_2ppn_SSI[trial], dnds_h_2ppn_SSD[trial],color = '0.5', alpha = 0.1, s = 10)
ax2.scatter(dnds_h_2ppn_SSI[0][2], dnds_h_2ppn_SSD[0][2],color = '0.5', alpha = 0.5, s = 30, label = r'$dN^h/dS^h$')
ax2.plot(line_x, line_x, linestyle = '--', color = 'k', marker = '' , label = 'x = y')
ax2.set_ylim(-0.05, 1.05)
ax2.set_xlim(-0.05, 1.05)
ax2.set_title('2PPN', fontsize = 16)

for x, y in zip(ppn_binned_SSI_dnds, ppn_binned_delta_dnds):
    ax5.errorbar(np.mean(x),np.mean(y), yerr = np.std(y), xerr = np.std(x), color = 'k', marker = 'o')
ax5.axhline(linestyle = '--', color = 'k')
ax5.set_yticks(np.arange(-0.2, 0.4, 0.1))
ax5.set_xlabel(r'expected $dN^h/dS^h$ under S-SI' , fontsize = 12)

## 1PEK  ##
for trial in range(50):
    ax3.scatter(dnds_h_1pek_SSI[trial], dnds_h_1pek_SSD[trial],color = '0.5', alpha = 0.1, s = 10)
ax3.scatter(dnds_h_1pek_SSI[0][2], dnds_h_1pek_SSD[0][2],color = '0.5', alpha = 0.5, s = 30, label = r'$dN^h/dS^h$')
ax3.plot(line_x, line_x, linestyle = '--', color = 'k', marker = '' , label = 'x = y')
ax3.set_ylim(-0.05, 1.05)
ax3.set_xlim(-0.05, 1.05)
ax3.set_title('1PEK', fontsize = 16)

for x, y in zip(pek_binned_SSI_dnds, pek_binned_delta_dnds):
    ax6.errorbar(np.mean(x),np.mean(y), yerr = np.std(y), xerr = np.std(x), color = 'k', marker = 'o')
ax6.axhline(linestyle = '--', color = 'k')
ax6.set_yticks(np.arange(-0.2, 0.4, 0.1))
ax6.set_xlabel(r'expected $dN^h/dS^h$ under S-SI' , fontsize = 12)

plt.savefig('../Figures/Figure7.png', dpi = 450, bbox_inches = 'tight')





