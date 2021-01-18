# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 14:56:53 2019

@author: libon
"""

import numpy as np
from scipy import stats


###########Calculate entropy 
def __shannon_entropy(c):
    p = c / np.sum(c)
    p = p[p > 0]
    return -np.sum(p * np.log2(p))
###########calculate the entropy


def __chaoShen(y):
    yx = y[y> 0]           # remove bins with zero counts
    n = sum(yx)             # total number of counts
    p = yx/n                # empirical frequencies
    f1 = sum(yx[yx ==1])    # number of singletons
    if f1 == n:
        f1 = n - 1
    C = 1 - f1/n            # estimated coverage
    gt = C*p                # coverage adjusted empirical frequencies (Good-Turing)
    denominator  = 1 - np.power((1-gt),n)
    H_cs = -np.sum(gt*np.log2(gt)/denominator)
    return H_cs

    
#######freedman_diaconis bin method
def freedman_diaconis(data):
    IQR  = stats.iqr(data,axis=0, rng=(25, 75), scale="raw", nan_policy="omit")
    N    = len(data)
    bw   = (2 * IQR) / np.power(N, 1/3)
    maxs, mins = np.max(data,axis = 0),np.min(data,axis = 0)
    datrng = maxs-mins
    result = np.array((datrng / bw) + 1).astype(int)
    return result #####amount of bin number for each column 
###################################


###############1D pdf estimation
def pdf_1D(Array_3D,i):
    Bin_num = freedman_diaconis(Array_3D)[i] 
    pdf_list = np.histogram(Array_3D[:,i],  bins = Bin_num)[0]
    zeros_1D = len(pdf_list[pdf_list != 0]) ##non-zero bins 
    
    H_mle = __shannon_entropy(pdf_list)
    H_cs = __chaoShen(pdf_list)
    H_mm = H_mle + (zeros_1D - 1)/(2*len(Array_3D))  
    return pdf_list, H_mle, H_cs, H_mm


#################2D pdf estimation 
def pdf_2D(Array,i,j):
    bin_s = freedman_diaconis(Array)
    pdf_list2d = np.histogram2d(Array[:,i], Array[:,j], 
                              bins = [bin_s[i], bin_s[j]])[0].flatten()
    zeros_2D = len(pdf_list2d[pdf_list2d != 0])
    
    H_mle2d = __shannon_entropy(pdf_list2d)
    H_cs2d = __chaoShen(pdf_list2d)
    H_mm2d = H_mle2d + (zeros_2D - 1)/(2*len(Array))
    return pdf_list2d, H_mle2d, H_cs2d, H_mm2d

###############3D pdf 
def pdf_3D_(Array,i,j,k):
    bin3d = freedman_diaconis(Array)
    pdf_list3d = np.histogramdd(Array[:,[i,j,k]], 
                              bins = [bin3d[i], bin3d[j], bin3d[k]])[0].flatten()
   
    zeros_3D = len(pdf_list3d[pdf_list3d != 0])
    H_mle3d = __shannon_entropy(pdf_list3d)
    H_cs3d = __chaoShen(pdf_list3d)
    H_mm3d = H_mle3d + (zeros_3D - 1)/(2*len(Array))
    return pdf_list3d, H_mle3d, H_cs3d, H_mm3d

###############4D pdf 
def pdf_4D_(Array,i,j,k,h):
    bin4d = freedman_diaconis(Array)
    pdf_list4d = np.histogramdd(Array[:,[i,j,k,h]], 
                              bins = [bin4d[i], bin4d[j], bin4d[k], bin4d[h]])[0].flatten()
    
    zeros_4D = len(pdf_list4d[pdf_list4d != 0])
    H_mle4d = __shannon_entropy(pdf_list4d)
    H_cs4d = __chaoShen(pdf_list4d)
    H_mm4d = H_mle4d + (zeros_4D - 1)/(2*len(Array))
    return pdf_list4d, H_mle4d, H_cs4d, H_mm4d










    

    

    




    




