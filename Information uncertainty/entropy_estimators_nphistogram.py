# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 14:56:53 2019

@author: libon
"""

import numpy as np
from scipy import stats
import copy 

###########Calculate entropy 
def shannonEntropy(c):
    c = copy.deepcopy(c)
    p = c / np.sum(c)
    p = p[p > 0]
    return -np.sum(p * np.log2(p))
###########calculate the entropy

    
#######freedman_diaconis bin method
def freedman_diaconis(data):
    data = copy.deepcopy(data)
    IQR  = stats.iqr(data,axis=0, rng=(25, 75), scale = 1, nan_policy="omit")
    N    = len(data)
    assert N == data.shape[0], "npts != data length, Stop and check!"
    bw   = (2 * IQR) / np.power(N, 1/3) ##binwidth = 2*IQR(x)/(n^1/3)
    maxs, mins = np.max(data,axis = 0),np.min(data,axis = 0)
    datrng = maxs-mins
    result = np.array((datrng / bw) + 1).astype(int) ##number of bins needed
    return result #####amount of bin number for each column 
##################################


###############1D pdf estimation
def pdf_1D(Array_3D,i):
    Array_3D = copy.deepcopy(Array_3D)
    assert len(Array_3D) == Array_3D.shape[0], "npts != data length, Stop and check!"
    Bin_num = freedman_diaconis(Array_3D)[i] 
    pdf_list = np.histogram(Array_3D[:,i],  bins = Bin_num)[0]
    zeros_1D = len(pdf_list[pdf_list != 0]) ##number of non-zero bins 
    H_mle = shannonEntropy(pdf_list) ##MLE entropy
    H_mm = H_mle + (zeros_1D - 1)/(2*len(Array_3D)) ##  acutually this correction
                                                    ##  may not as good as expected 
                                                    ##  users feel free to drop this correction
    return H_mm


#################2D pdf estimation 
def pdf_2D(Array,i,j):
    Array = copy.deepcopy(Array)
    assert len(Array) == Array.shape[0], "npts != data length, Stop and check!"
    bin_s = freedman_diaconis(Array)
    pdf_list2d = np.histogram2d(Array[:,i], Array[:,j], 
                              bins = [bin_s[i], bin_s[j]])[0].flatten()
    zeros_2D = len(pdf_list2d[pdf_list2d != 0])
    
    H_mle2d = shannonEntropy(pdf_list2d)
    H_mm2d = H_mle2d + (zeros_2D - 1)/(2*len(Array))
    return H_mm2d

###############3D pdf 
def pdf_3D(Array,i,j,k):
    Array = copy.deepcopy(Array)
    assert len(Array) == Array.shape[0], "npts != data length, Stop and check!"
    bin3d = freedman_diaconis(Array)
    pdf_list3d = np.histogramdd(Array[:,[i,j,k]], 
                              bins = [bin3d[i], bin3d[j], bin3d[k]])[0].flatten()
   
    zeros_3D = len(pdf_list3d[pdf_list3d != 0])
    H_mle3d = shannonEntropy(pdf_list3d)
    H_mm3d = H_mle3d + (zeros_3D - 1)/(2*len(Array))
    return H_mm3d

###############4D pdf 
def pdf_4D(Array,i,j,k,h):
    Array = copy.deepcopy(Array)
    assert len(Array) == Array.shape[0], "npts != data length, Stop and check!"
    bin4d = freedman_diaconis(Array)
    pdf_list4d = np.histogramdd(Array[:,[i,j,k,h]], 
                              bins = [bin4d[i], bin4d[j], bin4d[k], bin4d[h]])[0].flatten()
    
    zeros_4D = len(pdf_list4d[pdf_list4d != 0])
    H_mle4d = shannonEntropy(pdf_list4d)
    H_mm4d = H_mle4d + (zeros_4D - 1)/(2*len(Array))
    return H_mm4d










    

    

    




    




