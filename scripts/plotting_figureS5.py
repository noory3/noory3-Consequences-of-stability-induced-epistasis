#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 15:32:24 2020

@author: nooryoussef

Plotting bias and MSE
"""
import pandas as pd 
import seaborn as sns 
import matplotlib.pyplot as plt


df  = pd.read_csv("../Results/bias_and_MSE.csv")
sns.set_context("notebook")

f, axes = plt.subplots(1, 2, figsize = (8,4))
sns.boxplot(x="protein", y="bias", data=df, ax= axes[0], palette= ['0.9', '0.9', '0.9'], 
            order = ["1QHW", "2PPN", "1PEK"])
axes[0].set_ylim(-0.0001, 0.0004)
axes[0].set_xlabel("")
axes[0].set_ylabel("$bias^h$")

sns.boxplot(x="protein", y="MSE", data=df, ax= axes[1], palette= ['0.9', '0.9', '0.9'],
                        order = ["1QHW", "2PPN", "1PEK"])
axes[1].set_ylim(-0.00005, 0.0008)
axes[1].set_xlabel("")
axes[1].set_ylabel("$MSE^h$")

f.tight_layout(pad=2)
f.savefig("../Figures/bias_and_MSE.png")