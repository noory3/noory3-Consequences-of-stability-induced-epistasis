#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 13:13:18 2020

@author: nooryoussef

MLEs from M3 k =2 for all models and real proteins 
+ the number of amino acids per site 
"""
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import helpers

plt.rcParams["font.family"] = "Arial"
plt.rcParams["mathtext.fontset"] = 'cm'

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
    
    RealProp = Proportions(CountAAperSite(ReadSequences(path_to_alignments + protein +"_numbered.txt")))
    
    if protein == '2ppn':
        num_sites -= 1 
    elif protein == '1pek':
        num_sites = 279
        
    CSI = FullSimProp(path_to_alignments + protein + "_Ne2_C-SI")
    SSI = FullSimProp(path_to_alignments + protein + "_Ne2_S-SI")
    SSD = FullSimProp(path_to_alignments + protein + "_Ne2_S-SD")
    
    #mean number of sites with one, two, .., six amino acids averaged across trials
    one = []; two = []; three = []; four = []; five = []
    e_one = []; e_two = []; e_three = []; e_four = []; e_five = []
    for model in [CSI, SSI, SSD]:
        # mean # 
        one.append(np.mean(model, axis = 0)[0]*100)
        two.append(np.mean(model, axis = 0)[1]*100)
        three.append(np.mean(model, axis = 0)[2]*100)
        four.append(np.mean(model, axis = 0)[3]*100)
        five.append(sum(np.mean(model, axis = 0)[4:])*100)
    
        # std # 
        e_one.append(np.std(model, axis = 0)[0]*100)
        e_two.append(np.std(model, axis = 0)[1]*100)
        e_three.append(np.std(model, axis = 0)[2]*100)
        e_four.append(np.std(model, axis = 0)[3]*100)
        e_five.append(np.mean(np.std(model, axis = 0)[4:])*100)
    
    #add real PAP #     
    one.append(RealProp[0]*100)
    two.append(RealProp[1]*100)
    three.append(RealProp[2]*100)
    four.append(RealProp[3]*100)
    five.append(sum(RealProp[4:])*100)
    e_one.append(0.00)
    e_two.append(0.00)
    e_three.append(0.00)
    e_four.append(0.000)
    e_five.append(0.000)
    return(one, two, three, four, five, e_one, e_two, e_three, e_four, e_five)

#%%
RstPath = '../Results/Inference/'

#1qhw#
qhw_CSI   = pd.read_csv(RstPath + '1qhw_Ne2_C-SI_M3_k2.csv', sep = " ")
qhw_SSI   = pd.read_csv(RstPath + '1qhw_Ne2_S-SI_M3_k2.csv', sep = " ")
qhw_SSD   = pd.read_csv(RstPath + '1qhw_Ne2_S-SD_M3_k2.csv', sep = " ")

#2ppn#
ppn_CSI   = pd.read_csv(RstPath + '2ppn_Ne2_C-SI_M3_k2.csv', sep = " ")
ppn_SSI   = pd.read_csv(RstPath + '2ppn_Ne2_S-SI_M3_k2.csv', sep = " ")
ppn_SSD   = pd.read_csv(RstPath + '2ppn_Ne2_S-SD_M3_k2.csv', sep = " ")

#1pek#
pek_CSI   = pd.read_csv(RstPath + '1pek_Ne2_C-SI_M3_k2.csv', sep = " ")
pek_SSI   = pd.read_csv(RstPath + '1pek_Ne2_S-SI_M3_k2.csv', sep = " ") 
pek_SSD   = pd.read_csv(RstPath + '1pek_Ne2_S-SD_M3_k2.csv', sep = " ")

#%% create big data frame 
qhw_data = {'w' : list(qhw_CSI.w0) + list(qhw_CSI.w1) +  list(qhw_SSI.w0) +  list(qhw_SSI.w1) + list(qhw_SSD.w0) + list(qhw_SSD.w1) ,
          'p' : list(qhw_CSI.p0) + list(qhw_CSI.p1) +  list(qhw_SSI.p0) +  list(qhw_SSI.p1) + list(qhw_SSD.p0) + list(qhw_SSD.p1) ,
          'category' : ['1']*50+['2']*50 + ['1']*50+['2']*50 + ['1']*50+['2']*50 , 
          'model': ['C-SI']*100 + ['S-SI']*100 + ['S-SD']*100 }

ppn_data = {'w' : list(ppn_CSI.w0) + list(ppn_CSI.w1) +  list(ppn_SSI.w0) +  list(ppn_SSI.w1) + list(ppn_SSD.w0) + list(ppn_SSD.w1) ,
          'p' : list(ppn_CSI.p0) + list(ppn_CSI.p1) +  list(ppn_SSI.p0) +  list(ppn_SSI.p1) + list(ppn_SSD.p0) + list(ppn_SSD.p1) ,
          'category' : ['1']*50+['2']*50 + ['1']*50+['2']*50 + ['1']*50+['2']*50 , 
          'model': ['C-SI']*100 + ['S-SI']*100 + ['S-SD']*100 }

pek_data = {'w' : list(pek_CSI.w0) + list(pek_CSI.w1) +  list(pek_SSI.w0) +  list(pek_SSI.w1) + list(pek_SSD.w0) + list(pek_SSD.w1) ,
          'p' : list(pek_CSI.p0) + list(pek_CSI.p1) +  list(pek_SSI.p0) +  list(pek_SSI.p1) + list(pek_SSD.p0) + list(pek_SSD.p1) ,
          'category' : ['1']*50+['2']*50 + ['1']*50+['2']*50 + ['1']*50+['2']*50 , 
          'model': ['C-SI']*100 + ['S-SI']*100 + ['S-SD']*100 }

qhw_df = pd.DataFrame(qhw_data)
ppn_df = pd.DataFrame(ppn_data)
pek_df = pd.DataFrame(pek_data)

#%%#################################################################
#############          Plotting         ##########################
##################################################################
import seaborn as sns
colors = ["royalblue", "aliceblue", "red","mistyrose", "burlywood", "oldlace"]
line_colors = ["darkblue", "darkblue", "darkred", "darkred", "chocolate", "chocolate"]
    
def set_boxplot_colors(axis):
    for i,mybox in enumerate(axis.artists):
            mybox.set_facecolor(colors[i])
            mybox.set_edgecolor(line_colors[i])
            mybox.set_linewidth(1.5)
            for j in range(i*6,i*6+6):
                line = axis.lines[j]
                line.set_color(line_colors[i])
                line.set_mfc(line_colors[i])
                line.set_mec(line_colors[i])



barWidth = 0.1
names = ['C-SI', 'S-SI','S-SD', 'real' ]
x = [x for x in range(0, len(names))]

# Set position of bar on X axis
r1 = [-0.2, 0.8, 1.8, 2.8]
r2 = [x + barWidth for x in r1]
r3 = [x + barWidth for x in r2]
r4 = [x + barWidth for x in r3]
r5 = [x + barWidth for x in r4]

#%%#### Plott ####

fig, ((ax1, ax4, ax7), (ax2, ax5, ax8), (ax3, ax6, ax9))  = plt.subplots(3, 3, sharex = True, sharey = 'row', figsize = ((12,8)))
fig.subplots_adjust(wspace=0.15, hspace = 0.1)

### 1qhw ### 
df = qhw_df 
one, two, three, four, five, e_one, e_two, e_three, e_four, e_five = num_aa_per_site('1qhw')

#omega# 
sns.boxplot(x='model', y='w', hue='category', data=df, ax = ax1)
set_boxplot_colors(ax1)
ax1.set_ylabel(r'$\omega$', fontsize = 14)
ax1.set_xlabel('')
ax1.legend_.remove()
ax1.scatter([2.8,3.2], [0.007, 0.298], marker="o", facecolors= ['k', '0.9'], edgecolors="k")
ax1.set_title('1QHW \n tree length = 4.93', fontsize = 16)

#prop# 
sns.boxplot(x='model', y='p', hue='category', data=df, ax = ax2)
set_boxplot_colors(ax2)
ax2.set_ylabel(r'$p$', fontsize = 14)
ax2.set_xlabel('')
ax2.legend_.remove()
ax2.scatter([2.8,3.2], [0.715, 0.285], marker="o", facecolors= ['k', '0.9'], edgecolors="k")

#num aa per site # 
ax3.bar(r1, one,   yerr = e_one,   color=line_colors[::2] + ["k"], edgecolor= line_colors[::2] + ["k"], alpha = 0.2, width=barWidth)
ax3.bar(r2, two,   yerr = e_two,   color=line_colors[::2] + ["k"], edgecolor= line_colors[::2] + ["k"], alpha = 0.4, width=barWidth)
ax3.bar(r3, three, yerr = e_three, color=line_colors[::2] + ["k"], edgecolor= line_colors[::2] + ["k"], alpha = 0.6, width=barWidth)
ax3.bar(r4, four,  yerr = e_four,  color=line_colors[::2] + ["k"], edgecolor= line_colors[::2] + ["k"], alpha = 0.8, width=barWidth)
ax3.bar(r5, five,  yerr = e_five,  color=line_colors[::2] + ["k"], edgecolor= line_colors[::2] + ["k"], alpha = 1.0, width=barWidth)
ax3.set_ylim(0,85)
ax3.set_ylabel('% sites', fontsize =12)


### 2ppn ### 
df = ppn_df 
one, two, three, four, five, e_one, e_two, e_three, e_four, e_five = num_aa_per_site('2ppn')

#omega# 
sns.boxplot(x='model', y='w', hue='category', data=df, ax = ax4, fliersize=0)
set_boxplot_colors(ax4)
ax4.set_xlabel('')
ax4.scatter([2.8,3.2], [0.003, 0.086], marker="o", facecolors= ['k', '0.9'], edgecolors="k")
ax4.legend_.remove()
ax4.set_ylabel('')
ax4.set_title('2PPN\n tree length = 8.04', fontsize = 16)

#prop# 
sns.boxplot(x='model', y='p', hue='category', data=df, ax = ax5, fliersize=0)
set_boxplot_colors(ax5)
ax5.set_xlabel('')
ax5.legend_.remove()
ax5.scatter([2.8,3.2], [0.671, 0.329], marker="o", facecolors= ['k', '0.9'], edgecolors="k")
ax5.set_ylabel('')

#num aa per site # 
ax6.bar(r1, one,   yerr = e_one,   color=line_colors[::2] + ["k"], edgecolor= line_colors[::2] + ["k"], alpha = 0.2, width=barWidth)
ax6.bar(r2, two,   yerr = e_two,   color=line_colors[::2] + ["k"], edgecolor= line_colors[::2] + ["k"], alpha = 0.4, width=barWidth)
ax6.bar(r3, three, yerr = e_three, color=line_colors[::2] + ["k"], edgecolor= line_colors[::2] + ["k"], alpha = 0.6, width=barWidth)
ax6.bar(r4, four,  yerr = e_four,  color=line_colors[::2] + ["k"], edgecolor= line_colors[::2] + ["k"], alpha = 0.8, width=barWidth)
ax6.bar(r5, five,  yerr = e_five,  color=line_colors[::2] + ["k"], edgecolor= line_colors[::2] + ["k"], alpha = 1.0, width=barWidth)
ax6.set_ylim(0,85)
ax6.set_ylabel('', fontsize =12)

### 1pek ### 
df = pek_df 
one, two, three, four, five, e_one, e_two, e_three, e_four, e_five = num_aa_per_site('1pek')

#omega# 
sns.boxplot(x='model', y='w', hue='category', data=df, ax = ax7, fliersize=0)
set_boxplot_colors(ax7)
ax7.set_ylabel(r'$\omega$', fontsize = 12)
ax7.set_xlabel('')
ax7.scatter([2.8,3.2], [0.021, 0.235], marker="o", facecolors= ['k', '0.9'], edgecolors="k")
ax7.set_ylabel('')
ax7.set_title('1PEK\n tree length = 13.88', fontsize = 16)
custom_lines = [Line2D([0], [0], color='royalblue', marker='s', markeredgecolor= "darkblue",  linestyle='None', label = ""),
                Line2D([0], [0], color='aliceblue', marker='s', markeredgecolor= "darkblue",  linestyle='None', label = ""), 
                Line2D([0], [0], color='red'      , marker='s', markeredgecolor= "darkred",   linestyle='None', label = ""),
                Line2D([0], [0], color='mistyrose', marker='s', markeredgecolor= "darkred",   linestyle='None', label = ""), 
                Line2D([0], [0], color='burlywood', marker='s', markeredgecolor= "chocolate", linestyle='None', label = ""),
                Line2D([0], [0], color='oldlace'  , marker='s', markeredgecolor= "chocolate", linestyle='None', label = ""),
                Line2D([0], [0], color='black'    , marker='o', markeredgecolor= "k",         linestyle='None', label = r"$\omega_1$"), 
                Line2D([0], [0], color='0.9'      , marker='o', markeredgecolor= "k",         linestyle='None', label = r"$\omega_2$")]

ax7.legend(handles= custom_lines, loc='upper right',title="", ncol=4, columnspacing= 0.001, labelspacing = 0.1, handletextpad=0 , frameon = True, fancybox = True)

#prop# 
sns.boxplot(x='model', y='p', hue='category', data=df, ax = ax8, fliersize=0)
set_boxplot_colors(ax8)
ax8.set_ylabel('', fontsize = 12)
ax8.set_xlabel('')
ax8.scatter([2.8,3.2], [0.644, 0.356], marker="o", facecolors= ['k', '0.9'], edgecolors="k")
ax8.set_ylim(0,1.15)
custom_lines = [Line2D([0], [0], color='royalblue', marker='s', markeredgecolor= "darkblue",  linestyle='None', label = ""),
                Line2D([0], [0], color='aliceblue', marker='s', markeredgecolor= "darkblue",  linestyle='None', label = ""), 
                Line2D([0], [0], color='red'      , marker='s', markeredgecolor= "darkred",   linestyle='None', label = ""),
                Line2D([0], [0], color='mistyrose', marker='s', markeredgecolor= "darkred",   linestyle='None', label = ""), 
                Line2D([0], [0], color='burlywood', marker='s', markeredgecolor= "chocolate", linestyle='None', label = ""),
                Line2D([0], [0], color='oldlace'  , marker='s', markeredgecolor= "chocolate", linestyle='None', label = ""),
                Line2D([0], [0], color='black'    , marker='o', markeredgecolor= "k",         linestyle='None', label = r"$p_1$"), 
                Line2D([0], [0], color='0.9'      , marker='o', markeredgecolor= "k",         linestyle='None', label = r"$p_2$")]

ax8.legend(handles= custom_lines, loc='upper right',title="", ncol=4, columnspacing= 0.001, labelspacing = 0.1, handletextpad=0 , frameon = True, fancybox = True)

#num aa per site # 
ax9.bar(r1, one,   yerr = e_one,   color=line_colors[::2] + ["k"], edgecolor= line_colors[::2] + ["k"], alpha = 0.2, width=barWidth)
ax9.bar(r2, two,   yerr = e_two,   color=line_colors[::2] + ["k"], edgecolor= line_colors[::2] + ["k"], alpha = 0.4, width=barWidth)
ax9.bar(r3, three, yerr = e_three, color=line_colors[::2] + ["k"], edgecolor= line_colors[::2] + ["k"], alpha = 0.6, width=barWidth)
ax9.bar(r4, four,  yerr = e_four,  color=line_colors[::2] + ["k"], edgecolor= line_colors[::2] + ["k"], alpha = 0.8, width=barWidth)
ax9.bar(r5, five,  yerr = e_five,  color=line_colors[::2] + ["k"], edgecolor= line_colors[::2] + ["k"], alpha = 1.0, width=barWidth)
ax9.set_ylim(0,105)
ax9.set_ylabel('', fontsize =12)
custom_lines = [Line2D([0], [0], color='0.85', lw=6),
                Line2D([0], [0], color='0.70', lw=6),
                Line2D([0], [0], color='0.55', lw=6),
                Line2D([0], [0], color='0.40', lw=6),
                Line2D([0], [0], color='0.00', lw=6)]

ax9.legend(custom_lines, ["1", "2", "3", "4", r"$\geq$5"], loc='upper right',title="number of amino acids per site", ncol=3, labelspacing = 0.1, frameon = True, fancybox = True)

plt.xlim(-0.6,3.4)
plt.xticks([0,1,2,3], names)
[t.set_color(i) for (i,t) in zip(["darkblue", "darkred", "chocolate", "k"],ax3.xaxis.get_ticklabels())]
[t.set_color(i) for (i,t) in zip(["darkblue", "darkred", "chocolate", "k"],ax6.xaxis.get_ticklabels())]
[t.set_color(i) for (i,t) in zip(["darkblue", "darkred", "chocolate", "k"],ax9.xaxis.get_ticklabels())]



plt.savefig('../Figures/Figure4.png', dpi = 450, bbox_inches = 'tight')


