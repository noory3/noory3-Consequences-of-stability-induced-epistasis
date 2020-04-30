#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 14 08:34:20 2019

@author: N Youssef

Prints # significant trials and the average mles for models 
M0, M3 (k=2), M3 (k=3) and CLM3
"""
import pandas as pd 
import numpy as np

RstPath = '../Results/Inference/'

def omega1_less_than_omega2(df):
    '''
        Given the CLM3 results dataframe, enforces w1 < w2 
    '''
    for i, omega_1 in enumerate(df.w1):
        omega_2 = df.w2[i]
        if omega_1 > omega_2:
            df.w1[i] = omega_2
            df.w2[i] = omega_1
            df.p1[i] = 1 - df.p1[i]

    return(df)

def omega1_less_than_omega2_RaMoSS(df):
    '''
        Given the RaMoSS results dataframe, enforces w1 < w2 
    '''
    for i, omega_1 in enumerate(df.w1M3):
        omega_2 = df.w2M3[i]
        if omega_1 > omega_2:
            df.w1M3[i] = omega_2
            df.w2M3[i] = omega_1
            df.p1M3[i] = 1 - df.p1M3[i]

    return(df)
    
def BUSTED_rst(df):
    i=0
    w0_w1_w2_p0_p1 = np.zeros((50,5))
    for idx in range(len(df)):
        if df.iloc[idx].Model != "Unconstrained_model" and df.iloc[idx].Model != "Constrained_model":
            if str(df.iloc[idx+1].likelihood) == 'nan':   
                w0= df.iloc[idx+2].w0
                w1 = df.iloc[idx+2].w1
                w2 = df.iloc[idx+2].w2
                p0 = df.iloc[idx+2].p0
                p1 = df.iloc[idx+2].p1
                w0_w1_w2_p0_p1[i] = [w0, w1,w2,p0,p1]
                i+=1 
            elif str(df.iloc[idx+1].likelihood) != 'nan':
                w0= df.iloc[idx+1].w0
                w1 = df.iloc[idx+1].w1
                w2 = df.iloc[idx+1].w2
                p0 = df.iloc[idx+1].p0
                p1 = df.iloc[idx+1].p1
                w0_w1_w2_p0_p1[i] = [w0, w1,w2,p0,p1]
                i+=1 
    return(w0_w1_w2_p0_p1)
    
def BUSTED_sig(BUSTED):
    num_sig = sum(i < 0.005 for i in BUSTED['p-value'])
    rst = np.zeros((num_sig, 5))
    idx = [i for i,x in enumerate(BUSTED['p-value']) if x <= 1]
    r = 0
    for i in idx:
        x = BUSTED['p-value'][i]
        if x < 0.005:
            w0 =  BUSTED['w0'][i]
            w1 =  BUSTED['w1'][i]
            w2 =  BUSTED['w2'][i]
            p1 =  BUSTED['p0'][i]
            p2 =  BUSTED['p1'][i]
            rst[r] = [w0, w1, w2, p1, p2]
            r +=1
    return(rst)
 
#%%
protein = '1pek'
model = 'S-SI'
gen_model = protein + '_Ne2_' + model
print('\n\n######### gen model ' + str(gen_model) + '  #################')

M0        = pd.read_csv(RstPath + gen_model + '_M0.csv', sep = " ")
M3k2      = pd.read_csv(RstPath + gen_model + '_M3_k2.csv', sep = " ")
M3k3      = pd.read_csv(RstPath + gen_model + '_M3_k3.csv', sep = " ")
M3k4      = pd.read_csv(RstPath + gen_model + '_M3_k4.csv', sep = " ")
CLM3_null = pd.read_csv(RstPath + gen_model + '_CLM3_Null.csv', sep = " ")
CLM3_alt  = omega1_less_than_omega2(pd.read_csv(RstPath + gen_model + '_CLM3_Alt.csv', sep = " "))
BUSTED    = pd.read_csv(RstPath + gen_model + '_BUSTED_rst.csv', sep = " ")

# M0 # 
print ("M0")
print("w = " + "%.3f" % np.mean(M0.w))

#M3 k2#
k1vk2 = 2*(M3k2.LL - M0.LL)
print("\nM0 vs M3k2: " + str(sum(i > 5.99 for i in k1vk2)))

#Average mle's     
print("M3 k=2")
print("w1 = " + "%.3f" % np.mean(M3k2.w0[k1vk2>5.99]) + " (" + "%.3f" % np.mean(M3k2.w0) + ")")
print("w2 = " + "%.3f" % np.mean(M3k2.w1[k1vk2>5.99]) + " (" + "%.3f" % np.mean(M3k2.w1) + ")")
print("p1 = " + "%.3f" % np.mean(M3k2.p0[k1vk2>5.99]) + " (" + "%.3f" % np.mean(M3k2.p0) + ")")
 
#M3 k3#
k2vk3 = 2*(M3k3.LL - M3k2.LL)
print("\nM3k2 vs M3k3: " + str(sum(i > 5.99 for i in k2vk3)))

#Average mle's     
print("M3 k=3")
print("w1 = " + "%.3f" % np.mean(M3k3.w0[k2vk3>5.99]) + " (" + "%.3f" % np.mean(M3k3.w0)+ ")")
print("w2 = " + "%.3f" % np.mean(M3k3.w1[k2vk3>5.99]) + " (" + "%.3f" % np.mean(M3k3.w1)+ ")")
print("w3 = " + "%.3f" % np.mean(M3k3.w2[k2vk3>5.99]) + " (" + "%.3f" % np.mean(M3k3.w2)+ ")")
print("p1 = " + "%.3f" % np.mean(M3k3.p0[k2vk3>5.99]) + " (" + "%.3f" % np.mean(M3k3.p0)+ ")")
print("p2 = " + "%.3f" % np.mean(M3k3.p1[k2vk3>5.99]) + " (" + "%.3f" % np.mean(M3k3.p1)+ ")")

#M3 k4#
k3vk4 = 2*(M3k4.LL - M3k3.LL)
print("\nM3k3 vs M3k4: " + str(sum(i > 5.99 for i in k3vk4)))

#Average mle's     
print("M3 k=4")
print("w1 = " + "%.3f" % np.mean(M3k4.w0[k3vk4>5.99]) +  " (" + "%.3f" % np.mean(M3k4.w0)+ ")")
print("w2 = " + "%.3f" % np.mean(M3k4.w1[k3vk4>5.99]) +  " (" + "%.3f" % np.mean(M3k4.w1)+ ")")
print("w3 = " + "%.3f" % np.mean(M3k4.w2[k3vk4>5.99]) +  " (" + "%.3f" % np.mean(M3k4.w2)+ ")")
print("w4 = " + "%.3f" % np.mean(M3k4.w3[k3vk4>5.99]) +  " (" + "%.3f" % np.mean(M3k4.w3)+ ")")
print("p1 = " + "%.3f" % np.mean(M3k4.p0[k3vk4>5.99]) +  " (" + "%.3f" % np.mean(M3k4.p0)+ ")")
print("p2 = " + "%.3f" % np.mean(M3k4.p1[k3vk4>5.99]) +  " (" + "%.3f" % np.mean(M3k4.p1)+ ")")
print("p3 = " + "%.3f" % np.mean(M3k4.p2[k3vk4>5.99]) +  " (" + "%.3f" % np.mean(M3k4.p2)+ ")")

#CLM3#
M3k2vCLM3 =  2*(CLM3_null.LL - CLM3_alt.LL)
print("\nM3k2 vs CLM3: " + str(sum(i > 2.71 for i in M3k2vCLM3)))

#Average mle's     
print("CLM3")
print("w1 = "    + "%.3f" % np.mean(CLM3_alt.w1[M3k2vCLM3>2.71]) +  " (" + "%.3f" % np.mean(CLM3_alt.w1)+ ")")
print("w2 = "    + "%.3f" % np.mean(CLM3_alt.w2[M3k2vCLM3>2.71]) +  " (" + "%.3f" % np.mean(CLM3_alt.w2)+ ")")
print("p1 = "    + "%.3f" % np.mean(CLM3_alt.p1[M3k2vCLM3>2.71]) +  " (" + "%.3f" % np.mean(CLM3_alt.p1)+ ")")
print("delta = " + "%.3f" % np.mean(CLM3_alt.delta[M3k2vCLM3>2.71]) +  " (" + "%.3f" % np.mean(CLM3_alt.delta)+ ")")

#BUSTED# 
print("\nBUSTED: " + str(sum(i < 0.005 for i in BUSTED['p-value'])))
#    BUSTED = BUSTED_rst(BUSTED)
if sum(i < 0.005 for i in BUSTED['p-value']) > 0:
    BUSTED = BUSTED_sig(BUSTED)
    print("BUSTED")
    print("w1 = " + "%.3f" % np.mean(BUSTED[:,0]))
    print("w2 = " + "%.3f" % np.mean(BUSTED[:,1]))
    print("w3 = " + "%.3f" % np.mean(BUSTED[:,2]))
    print("p1 = " + "%.3f" % np.mean(BUSTED[:,3]))
    print("p2 = " + "%.3f" % np.mean(BUSTED[:,4]))
    

