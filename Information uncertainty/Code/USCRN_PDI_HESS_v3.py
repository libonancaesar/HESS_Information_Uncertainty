# -*- coding: utf-8 -*-
"""
Created on Sat Jul 11 19:08:10 2020

@author: libon
"""
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import concurrent.futures
# import Histogram_Entropy_Estimator as pdf
import entropy_estimators_nphistogram as pdf
from functools import partial
from scipy.stats import pearsonr,spearmanr, ttest_ind, ttest_rel
from dateutil.parser import parse
from datetime import datetime
from datetime import timedelta
from sklearn.metrics import mean_squared_error

Box = 'C:/Users/libon/Box/' ###box folder 
UScrn_sm = Box + 'USCRN_SM/' ###USCRN soil moisture folder 
plt.rcParams["figure.dpi"] = 300
plt.style.use('seaborn')
def get_Cleaned(product): 
    UTC_T = datetime(2000,1,1,0,0)
    if product == '9km':
        path = Box + 'USCRM_SMAP_sm/SMAP_crn9_vwc_temp.csv'
        names = ['site','smdca9km', 'tbh9km','tbv9km','temp9km','vwc9km', 'time', 'DOY9']
    elif product == '36km':
        path = Box + 'USCRM_SMAP_sm/SMAP_crn36_vwc_temp.csv'
        names = ['site','smdca36km', 'tbh36km','tbv36km','temp36km','vwc36km','lc_Int','lc_frac',
                 'lc_class','time', 'DOY36']
    SMAP = pd.read_csv(path, index_col= 0)
    
    #############convert SMAP time to UTC time to match with USCRN observations
    SMAP['utc'] = [UTC_T + timedelta(seconds = int(k)) for k in SMAP['utc']]
    ################
    SMAP['time'] = SMAP['utc'].dt.round('H')
    SMAP['DOY'] = SMAP['time'].dt.dayofyear
    SMAP_cleaned = SMAP.drop(['date','utc'],axis = 1)
    SMAP_cleaned.columns = names
    sites =  SMAP_cleaned['site'].unique()
    return SMAP_cleaned, sites

def NSE(sim, insitu):
    return 1 - np.sum((insitu - sim)**2)/np.sum((insitu - np.mean(insitu))**2) 



def getInformation(smap36km, smap9km, sitename):
        
    smap9km = smap9km[smap9km['site'] == sitename]
    smap36km = smap36km[smap36km['site'] == sitename]
    vwc9km = [smap9km[smap9km['DOY9'] == i]['vwc9km'].values[0] for i in smap9km['DOY9'].unique()]
    vwc36km = [smap36km[smap36km['DOY36'] == i]['vwc36km'].values[0] for i in smap36km['DOY36'].unique()]
    
    datahub = smap9km.merge(smap36km, on = 'time').drop(['site_x', 'site_y', 'DOY9', 'DOY36'], axis = 1)
    landcover = datahub['lc_class'].unique()[0]
    # return datahub
    ##USCRN SM
    data = pd.read_csv(UScrn_sm + sitename + '.txt',sep='\s+',header=None) 
    lat = data.iloc[0,3]
    long = data.iloc[0,2]
    
    data = data.loc[:,[1,4]]
    data.columns = ['time','sm5cm']
    data['time'] =  pd.to_datetime(data['time'], format="%Y%m%d%H")
    data = data[~data.sm5cm.isin([-9,-99,-999,-9999])]
     

    ##col0:'smdca9km', col1: 'tbh9km', col2:'tbv9km', col3:'temp9km'
    ##col4:'smdca36km', col5:'tbh36km', col6: 'tbv36km', col7: 'temp36km'
    ##col8:'sm5cm' (in situ)
    
    

    smap_crnSM = pd.merge(datahub,data, on = 'time').drop(['time','vwc9km','vwc36km',
                                                           'lc_Int','lc_frac','lc_class'],axis=1)
    
    print(smap_crnSM.columns.tolist(), landcover)
    if len(smap_crnSM) >= 500:
        # plt.figure()
        plt.scatter(smap_crnSM.smdca36km, smap_crnSM.smdca9km)
        # plt.scatter(smap_crnSM.sm5cm, smap_crnSM.smdca9km, label = '9')
        
        plt.xlabel('DCA 36KM')
        plt.ylabel('DCA 9KM')
        plt.plot([0, 0.45], [0, 0.45])
        print(pearsonr(smap_crnSM.smdca36km, smap_crnSM.smdca9km)[0])
        # print(smap_crnSM.head())
        Mean_vwc9km = sum(vwc9km)/len(vwc9km)
        Mean_vwc36km = sum(vwc36km)/len(vwc36km)
        N = np.log2(len(smap_crnSM))
        pear_R36 = pearsonr(smap_crnSM['smdca36km'],smap_crnSM['sm5cm'])[0]
        pear_R9 = pearsonr(smap_crnSM['smdca9km'],smap_crnSM['sm5cm'])[0]
        RMSE_36 = np.sqrt(mean_squared_error(smap_crnSM['smdca36km'],smap_crnSM['sm5cm']))
        RMSE_9 = np.sqrt(mean_squared_error(smap_crnSM['smdca9km'],smap_crnSM['sm5cm']))
        
        
        
        ########1D Entropy 
        H_dca_9 = pdf.pdf_1D(smap_crnSM.values,0)[-1]/N ## dca 9 km
        H_h_9 = pdf.pdf_1D(smap_crnSM.values,1)[-1]/N## h 9km
        H_v_9 = pdf.pdf_1D(smap_crnSM.values,2)[-1]/N ## v 9km
        
        H_dca_36 =  pdf.pdf_1D(smap_crnSM.values, 4)[-1]/N ## dca 36 km 
        H_h_36 = pdf.pdf_1D(smap_crnSM.values, 5)[-1]/N ##h 36 km 
        H_v_36 = pdf.pdf_1D(smap_crnSM.values, 6)[-1]/N ##v 36 km
        H_crn =  pdf.pdf_1D(smap_crnSM.values, 8)[-1]/N ## UScrn in situ  
        
        H_temp36 = pdf.pdf_1D(smap_crnSM.values, 7)[-1]/N
        H_temp9 = pdf.pdf_1D(smap_crnSM.values, 3)[-1]/N
        
        ########2D Entropy 
        H_hv_9 =  pdf.pdf_2D(smap_crnSM.values, 1, 2)[-1]/N
        H_h_dca_9 = pdf.pdf_2D(smap_crnSM.values, 0, 1)[-1]/N
        H_v_dca_9 = pdf.pdf_2D(smap_crnSM.values, 0, 2)[-1]/N
        H_dca_crn_9 = pdf.pdf_2D(smap_crnSM.values, 0, 8)[-1]/N
       
        H_hv_36 =  pdf.pdf_2D(smap_crnSM.values, 5, 6)[-1]/N
        H_h_dca_36 = pdf.pdf_2D(smap_crnSM.values, 4, 5)[-1]/N
        H_v_dca_36 = pdf.pdf_2D(smap_crnSM.values, 4, 6)[-1]/N
        H_dca_crn_36 = pdf.pdf_2D(smap_crnSM.values, 4 ,8)[-1]/N
        
        H_temp_crn_36 = pdf.pdf_2D(smap_crnSM.values, 7,8)[-1]/N
        H_temp_crn_9 = pdf.pdf_2D(smap_crnSM.values, 3,8)[-1]/N
        
        H_temp_dca_36 = pdf.pdf_2D(smap_crnSM.values, 4, 7)[-1]/N
        H_temp_dca_9 = pdf.pdf_2D(smap_crnSM.values, 0, 3)[-1]/N
        
        #######3D Entropy 
        H_hv_dca_9 = pdf.pdf_3D_(smap_crnSM.values, 0, 1, 2)[-1]/N    
        H_hv_crn_9 =  pdf.pdf_3D_(smap_crnSM.values, 1, 2, 8)[-1]/N
        
        H_hv_dca_36 =  pdf.pdf_3D_(smap_crnSM.values, 4, 5, 6)[-1]/N
        H_hv_crn_36 =  pdf.pdf_3D_(smap_crnSM.values, 5, 6, 8)[-1]/N
        
        H_hv_temp9 = pdf.pdf_3D_(smap_crnSM.values, 1, 2, 3)[-1]/N
        H_hv_temp36 = pdf.pdf_3D_(smap_crnSM.values, 5, 6, 7)[-1]/N
        

        
        #######4D Entropy
        H_hv_temp_crn9 = pdf.pdf_4D_(smap_crnSM.values, 1, 2, 3, 8)[-1]/N
        H_hv_temp_crn36 = pdf.pdf_4D_(smap_crnSM.values, 5, 6, 7, 8)[-1]/N
        
        H_hv_temp_dca9 = pdf.pdf_4D_(smap_crnSM.values, 0, 1, 2, 3)[-1]/N
        H_hv_temp_dca36 = pdf.pdf_4D_(smap_crnSM.values, 4, 5, 6, 7)[-1]/N
        
        
            ##################conditional mutual information 
        I_hvSM_given_Temp36 = H_hv_temp36 + H_temp_crn_36 - H_hv_temp_crn36 - H_temp36
        I_hvSM_given_Temp9 = H_hv_temp9 + H_temp_crn_9 - H_hv_temp_crn9 - H_temp9
        
        I_hvdca_given_Temp9 = H_hv_temp9 + H_temp_dca_9 - H_hv_temp_dca9 - H_temp9
        I_hvdca_given_Temp36 = H_hv_temp36 + H_temp_dca_36 - H_hv_temp_dca36 - H_temp36
        
        I_dca_givenhv36 = H_hv_temp36 + H_hv_dca_36 - H_hv_temp_dca36 - H_hv_36
        I_dca_givenhv9 = H_hv_temp9 + H_hv_dca_9 - H_hv_temp_dca9 - H_hv_9
        ####Mutual Information 3 Variable
        I_hv_dca_9 = H_hv_9 + H_dca_9 - H_hv_dca_9 
        I_hv_crn_9 =  H_hv_9 + H_crn - H_hv_crn_9
        
        I_hv_dca_36 = H_hv_36 + H_dca_36 - H_hv_dca_36 
        I_hv_crn_36 =  H_hv_36 + H_crn - H_hv_crn_36
        ###Mutual Information 2 Variable
        I_hv_9 = H_h_9 + H_v_9 - H_hv_9
        I_h_dca_9 = H_h_9 + H_dca_9 - H_h_dca_9
        I_v_dca_9 = H_v_9 + H_dca_9 - H_v_dca_9
        I_crn_dca_9 = H_crn + H_dca_9 - H_dca_crn_9 
        
        I_dca_temp36 = H_dca_36 + H_temp36 - H_temp_dca_36
        I_dca_temp9 = H_dca_9 + H_temp9 - H_temp_dca_9
        
        I_hv_36 = H_h_36 + H_v_36 - H_hv_36
        I_h_dca_36 = H_h_36 + H_dca_36 - H_h_dca_36
        I_v_dca_36 = H_v_36 + H_dca_36 - H_v_dca_36
        I_crn_dca_36 = H_crn + H_dca_36 - H_dca_crn_36 
        
        ##Mutual Information 4 varibles 
        I_hv_temp_crn9 = H_crn + H_hv_temp9 - H_hv_temp_crn9
        I_hv_temp_crn36 = H_crn + H_hv_temp36 - H_hv_temp_crn36
        
        I_hv_temp_dca9 = H_dca_9 + H_hv_temp9 - H_hv_temp_dca9
        I_hv_temp_dca36 = H_dca_36 + H_hv_temp36 - H_hv_temp_dca36
        
    
        
        
        #####PID
        II_9 = H_hv_9 + H_h_dca_9 + H_v_dca_9 - (H_h_9 + H_v_9 + H_dca_9) - H_hv_dca_9
        R_MMI_9 = min(I_h_dca_9,I_v_dca_9)
        R_min_9 = max(0,-II_9)
        I_s_9 = I_hv_9/min(H_h_9,H_v_9)
        R_s_9 = R_min_9 + I_s_9*(R_MMI_9-R_min_9)
        
        U_h_9 = I_h_dca_9 - R_s_9
        U_v_9 = I_v_dca_9 - R_s_9
        S_9 = II_9 + R_s_9
        
        II_36 = H_hv_36 + H_h_dca_36 + H_v_dca_36 - (H_h_36 + H_v_36 + H_dca_36) - H_hv_dca_36
        R_MMI_36 = min(I_h_dca_36,I_v_dca_36)
        R_min_36 = max(0,-II_36)
        I_s_36 = I_hv_36/min(H_h_36,H_v_36)
        R_s_36 = R_min_36 + I_s_36*(R_MMI_36-R_min_36)
        
        U_h_36 = I_h_dca_36 - R_s_36
        U_v_36 = I_v_dca_36 - R_s_36
        S_36 = II_36 + R_s_36
        
        # ########vectorized 
        # II_36 = H_hv_36 + H_h_dca_36 + H_v_dca_36 - (H_h_36 + H_v_36 + H_dca_36) - H_hv_dca_36
        # R_MMI_36 = min(I_h_dca_36,I_v_dca_36)
        # R_min_36 = max(0,-II_36)
        # I_s_36 = I_hv_36/min(H_h_36,H_v_36)
        # R_s_36 = R_min_36 + I_s_36*(R_MMI_36-R_min_36)
        
        # U_h_36 = I_h_dca_36 - R_s_36
        # U_v_36 = I_v_dca_36 - R_s_36
        # S_36 = II_36 + R_s_36         
        print('finished working on:',sitename)
        print('------------------------------------')
        return [sitename, lat, long, H_crn, H_h_9, H_v_9, H_temp9, H_dca_9, 
                H_hv_9, I_crn_dca_9, I_dca_temp9, I_hv_crn_9, I_hv_dca_9, I_hv_temp_crn9,
                I_hv_temp_dca9,I_hvSM_given_Temp9, I_hvdca_given_Temp9, I_dca_givenhv9,
                U_h_9, U_v_9, S_9, R_s_9, pear_R9, RMSE_9, Mean_vwc9km, 
                H_h_36, H_v_36,H_temp36, H_dca_36, H_hv_36, 
                I_crn_dca_36, I_dca_temp36, I_hv_crn_36, I_hv_dca_36,I_hv_temp_crn36,I_hv_temp_dca36,
                I_hvSM_given_Temp36, I_hvdca_given_Temp36, I_dca_givenhv36,
                U_h_36, U_v_36, S_36, R_s_36, pear_R36, RMSE_36,
                Mean_vwc36km,landcover, len(smap_crnSM)]
        
smap9, site9 = get_Cleaned('9km')
smap36, site36 = get_Cleaned('36km')
common_sets = [i for i in set.intersection(set(site9),set(site36))]

Info_Metric = [getInformation(smap36, smap9, i) for i in common_sets]
 
Entropy_MI = pd.DataFrame([i for i in Info_Metric if i is not None], 
                          columns = ['sitename','lat','long','H_in_situ','H_h_9','H_v_9','H_temp9','H_dca9', 
                                     'H_hv_9','I_crn_dca9', 'I_dca_temp9','I_hv_crn9', 'I_hv_dca_9','I_hv_temp_crn9',
                                     'I_hv_temp_dca9', 'I_hvSM_given_Temp9', 'I_hvdca_given_Temp9','I_dca_givenhv9',
                                      'U_h_9', 'U_v_9','S_9', 'R_9', 'Pear_R9','RMSE_9','vwc9km', 
                                      'H_h_36','H_v_36','H_temp36','H_dca36', 'H_hv_36',
                                      'I_crn_dca36', 'I_dca_temp36','I_hv_crn36', 'I_hv_dca_36', 'I_hv_temp_crn36','I_hv_temp_dca36',
                                      'I_hvSM_given_Temp36','I_hvdca_given_Temp36', 'I_dca_givenhv36',
                                      'U_h_36', 'U_v_36','S_36', 'R_36', 'Pear_R36', 'RMSE_36', 
                                      'vwc36km','landcover','len'])

Entropy_MI.to_csv(Box + 'USCRM_SMAP_sm/HESS_PID_MM.csv')
data = pd.read_csv(Box + 'USCRM_SMAP_sm/HESS_PID_MM.csv')



    




















