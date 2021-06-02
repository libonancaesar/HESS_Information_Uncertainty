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



##########generate 1d tuplelist 
def Grid_space(array1d,Bin_num):
    Grid = np.linspace(min(array1d),
                       max(array1d),
                       num = Bin_num+1, 
                       endpoint=True)
    Grid[-1] = Grid[-1]*1.1
    tuple_lists = [[Grid[i],Grid[i+1]]
                    for i in np.arange(len(Grid)) 
                    if i<=(len(Grid)-2)]
    return tuple_lists
#########generate 1d tuple list 
    
#######freedman_diaconis bin method
def freedman_diaconis(data):
    IQR  = stats.iqr(data,axis=0, rng=(25, 75), scale="raw", nan_policy="omit")
    N    = len(data)
    bw   = (2 * IQR) / np.power(N, 1/3)
    maxs, mins = np.max(data,axis = 0),np.min(data,axis = 0)
    datrng = maxs-mins
    result = np.array((datrng / bw) + 1).astype(int)
    # print(result)
    return result #####amount of bin number for each column 
###################################
   
    
###############1D pdf estimation
def pdf_1D(Array_3D,i):
    Bin_num = freedman_diaconis(Array_3D)[i] 
    tuple_lists = Grid_space(Array_3D[:,i],Bin_num)
    pdf_list = [0]*Bin_num
    for k in np.arange(len(tuple_lists)):
        for j in Array_3D[:,i]:
            if (tuple_lists[k][0] <= j < tuple_lists[k][1]):
                pdf_list[k] = pdf_list[k]+1
    m = len([i for i in pdf_list if i!=0])
    return pdf_list, __shannon_entropy(pdf_list) + (m-1)/(2*len(Array_3D))
###############1D pdf estimation 


#################2D pdf estimation 
def pdf_2D(Array,i,j):
    bin_i = freedman_diaconis(Array)[i]
    bin_j = freedman_diaconis(Array)[j]
    ######### first dimension 
    tuple_list_i = Grid_space(Array[:,i], bin_i)
    tuple_list_j = Grid_space(Array[:,j], bin_j)
    #############
    import itertools
    list_2d = list(itertools.product(tuple_list_i,tuple_list_j))
    pdf_2d_ = [0]*len(list_2d)
    for index in np.arange(len(list_2d)):
        for data in Array:
            if (
                    list_2d[index][0][0] <= data[i] < list_2d[index][0][1] and 
                    list_2d[index][1][0] <= data[j] < list_2d[index][1][1]
                    
                ):
                pdf_2d_[index] = pdf_2d_[index]+1
    m2 = len([i for i in pdf_2d_ if i!=0])
    return pdf_2d_, __shannon_entropy(pdf_2d_) + (m2-1)/(2*len(Array))
#######################################

###############3D pdf 
def pdf_3D_(Array,i,j,k):
    bins = freedman_diaconis(Array)
    t1 = Grid_space(Array[:,i], bins[i])
    t2 = Grid_space(Array[:,j], bins[j])
    t3 = Grid_space(Array[:,k], bins[k])
    import itertools
    list_3d = list(itertools.product(t1,t2,t3))
    pdf_3d__ = [0]*len(list_3d)
    
    for index in np.arange(len(list_3d)):
        for data in Array:
            if (
                    list_3d[index][0][0] <= data[i] < list_3d[index][0][1] and 
                    list_3d[index][1][0] <= data[j] < list_3d[index][1][1] and 
                    list_3d[index][2][0] <= data[k] < list_3d[index][2][1]
                ):
                pdf_3d__[index] = pdf_3d__[index]+1
    m3 = len([i for i in pdf_3d__ if i!=0])
    return pdf_3d__, __shannon_entropy(pdf_3d__)+ (m3-1)/(2*len(Array))
#######################################
    

###############4D pdf 
def pdf_4D_(Array,i,j,k,h):
    bins = freedman_diaconis(Array)
    t1 = Grid_space(Array[:,i], bins[i])
    t2 = Grid_space(Array[:,j], bins[j])
    t3 = Grid_space(Array[:,k], bins[k])
    t4 = Grid_space(Array[:,h], bins[h])
    import itertools
    list_4d = list(itertools.product(t1,t2,t3,t4))
    pdf_4d__ = [0]*len(list_4d)
    
    for index in np.arange(len(list_4d)):
        for data in Array:
            if (
                    list_4d[index][0][0] <= data[i] < list_4d[index][0][1] and 
                    list_4d[index][1][0] <= data[j] < list_4d[index][1][1] and 
                    list_4d[index][2][0] <= data[k] < list_4d[index][2][1] and 
                    list_4d[index][3][0] <= data[h] < list_4d[index][3][1]
                    
                ):
                pdf_4d__[index] = pdf_4d__[index]+1
    m4 = len([i for i in pdf_4d__ if i!=0])
    return pdf_4d__, __shannon_entropy(pdf_4d__) + (m4-1)/(2*len(Array))
#######################################










    

    

    




    




