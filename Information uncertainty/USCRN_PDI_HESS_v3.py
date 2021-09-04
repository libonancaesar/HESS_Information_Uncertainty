# -*- coding: utf-8 -*-
"""
Created on Sat Jul 11 19:08:10 2020

@author: libon
"""
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import entropy_estimators_nphistogram as pdf
from texttable import Texttable
from scipy.stats import pearsonr
from datetime import datetime
from datetime import timedelta

# from sklearn.metrics import mean_squared_error

Box = 'C:/Users/libon/Box/'
UScrn_sm = Box + 'USCRN_SM/'
colNames = ['site','smdca36km', 
            'tbh36km', 'tbv36km',
            'temp36km', 'vwc36km',
            'lc_Int','lc_frac',
            'lc_class','time', 'DOY36']
####processing the SMAP datasets Set this to 2000 "Noon" UTC time
####We explicitly say that the soil moisture retrieved at the tb time
####refer to 
####https://nsidc.org/sites/nsidc.org/files/technical-references/PSD_SPL2SMP_V7.pdf
###4.6.6 tb_time_seconds

UTC = datetime(2000,1,1,12,00) 

SMAP36km = pd.read_csv(Box + 'USCRM_SMAP_sm/SMAP_36km_loop_2021-05-27.csv', 
                   index_col= 0)

SMAP36km['utc'] = [UTC + timedelta(seconds = int(k)) for k in SMAP36km['utc']] 

###around the UTC time to the nearest hour 
SMAP36km['time'] = SMAP36km['utc'].dt.round('H')
SMAP36km['DOY'] = SMAP36km['time'].dt.dayofyear
SMAP_cleaned = SMAP36km.drop(['date','utc'],axis = 1)
SMAP_cleaned.columns = colNames
sites =  SMAP_cleaned['site'].unique()
lcIntUnique = SMAP_cleaned['lc_Int'].unique()
lcUnique = SMAP_cleaned['lc_class'].unique()

metrics = []

for i in sites:
    ###Load SMAP data product 
    print("Working on {}".format(i))
    SMAP_dt = SMAP_cleaned.loc[SMAP_cleaned['site'] == i, 
                              ['time','smdca36km', 
                               'tbh36km', 'tbv36km', 
                               'temp36km', 'lc_class']] 
    
    siteLandCover = list(SMAP_dt['lc_class'].unique())
 
    SMAP_dt = SMAP_dt.drop(columns = ['lc_class'])
    
    assert len(siteLandCover) == 1, "{} has multiple landcover stop and check data".format(i)
    
    sLC = siteLandCover[0]
    
    CRN = pd.read_csv(UScrn_sm + i + '.txt',sep='\s+',header=None) 
    ### The col 2 and col 3 are the site Lon and Lat
    ### col 1 and col 4 are the UTC time and 5cm soil moisture
    assert len(CRN.iloc[:,2].unique()) == 1
    assert len(CRN.iloc[:,3].unique()) == 1
    
    siteLat = CRN.iloc[0,3] 
    siteLon =  CRN.iloc[0,2]
    assert siteLon < 0, "site Long col # is 2, then SM 5cm col # is 4!"
    assert siteLat > 0, "site Long col # is 3, then SM 5cm col # is 4!"
    insitu = CRN.loc[:,[1,4]]

    insitu = insitu.rename(columns = {1:'time',4:'sm5cm'})
    ###In situ soil moisture are indeed in UTC times
    insitu['time'] =  pd.to_datetime(insitu['time'], format="%Y%m%d%H")
    insitu = insitu.loc[~insitu['sm5cm'].isin([-9,-99,-999,-9999])]
    
    print(insitu.drop_duplicates().shape[0], insitu.shape[0])
    print(SMAP_dt.drop_duplicates().shape[0], SMAP_dt.shape[0])

    SMAP_CRN =  pd.merge(insitu, SMAP_dt, on = 'time', how= 'inner')
    SMAP_CRN.drop_duplicates(subset=None, keep='first', 
                          inplace=True, ignore_index=True)
    SMAP_CRN = SMAP_CRN.dropna()
    
    minTime = min(SMAP_CRN['time'])
    maxTime = max(SMAP_CRN['time'])
    print(maxTime)
    SMAP_CRN = SMAP_CRN.drop(columns= ['time'])
    
    SMAP_CRN = SMAP_CRN[(SMAP_CRN['smdca36km'] > 0) & (SMAP_CRN['sm5cm'] > 0)]
    if False:
        plt.figure()
        plt.plot(SMAP_CRN['smdca36km'],label = 'SMAP36KM')
        plt.plot(SMAP_CRN['sm5cm'],label = 'CRN in situ')
        plt.legend()
        plt.title(i)
    
    numPoints = SMAP_CRN.shape[0]
    '''double check if the site can be dominated by different landcovers
        after checking that one site only dominated by one landcover!!!
    '''
    ###Load the LAI dataset 
    LAI = pd.read_csv("C:/Users/libon/Box/USCRM_SMAP_sm/USCRN_36KM_LAI/" +
                      i + '.csv')
    LAI = LAI[['date', 'range', 'mean', 'standard_deviation', 'pixels_pass_qa',
                'per_cent_pixels_pass_qa']]
    LAI = LAI.dropna()
    goodQualityLAI = LAI[(LAI['per_cent_pixels_pass_qa'] >= 70)].copy()
    goodQualityLAI['datetime'] = pd.to_datetime(goodQualityLAI['date'], format = '%Y-%m-%d')
    goodQualityLAI = goodQualityLAI[(goodQualityLAI['datetime'] >= minTime) &
                                    (goodQualityLAI['datetime'] <= maxTime)]
    ####standard deviation and mean of the mean LAI 
    std_LAI = goodQualityLAI['standard_deviation'].mean()
    mean_LAI = goodQualityLAI['mean'].mean()
   
    N = np.log2(SMAP_CRN.shape[0])     
    ####getting the correlation coefficents 
    corrSoilMoisture = pearsonr(SMAP_CRN['smdca36km'],SMAP_CRN['sm5cm'])[0] ##pearson R(in situ, DCA)
    corrHandV = pearsonr(SMAP_CRN['tbv36km'],SMAP_CRN['tbh36km'])[0]        ##pearson R(tbh,     tbv)
    rangeHV = np.mean(np.abs(SMAP_CRN['tbh36km'] - SMAP_CRN['tbv36km']))
    
    '''tbv is "GREATER" than tbh  
        if we dont remember which is which 
                  Data fields of 
                "SMAP_CRN.values"
      *===================================*
      | col0:'insitu'   | col1: 'dca36km', 
        -----------------------------------*
      | col2:'tbh36'    | col3:'tbv36'
        -----------------------------------*
      | col4:'temp36km' |
      *===================================*
    '''    

##################################

    
    H_insitu =  pdf.pdf_1D(SMAP_CRN.to_numpy(copy= True), 0)/N ## Entropy of H(in situ) 



    H_DCA =  pdf.pdf_1D(SMAP_CRN.to_numpy(copy= True), 1)/N    ## Entropy of H(dca)  36 km 
    H_Tbh = pdf.pdf_1D(SMAP_CRN.to_numpy(copy= True), 2)/N     ## Entropy of H(Tbh)  36 km 
    H_Tbv = pdf.pdf_1D(SMAP_CRN.to_numpy(copy= True), 3)/N     ## Entropy of H(Tbv) 36 km
    H_Teff = pdf.pdf_1D(SMAP_CRN.to_numpy(copy= True), 4)/N    ## Entropy of H(Teff) 36 km
    
    ##bivariate entropy
    H_h_v =  pdf.pdf_2D(SMAP_CRN.to_numpy(copy= True), 2, 3)/N        ## Joint entropy of H(tbh,tbv)
    H_h_dca = pdf.pdf_2D(SMAP_CRN.to_numpy(copy= True), 1, 2)/N       ## Joint entropy of H(tbh,DCA)
    H_v_dca = pdf.pdf_2D(SMAP_CRN.to_numpy(copy= True), 1, 3)/N       ## Joint entropy of H(tbv,DCA)
    H_dca_insitu = pdf.pdf_2D(SMAP_CRN.to_numpy(copy= True), 0, 1)/N  ## Joint entropy of H(DCA,in situ)
  
    ##Trivariate entropy
    H_h_v_dca =  pdf.pdf_3D(SMAP_CRN.to_numpy(copy= True), 1, 2, 3)/N    ## Joint entropy of H(tbh,tbv,DCA)
    H_h_v_temp = pdf.pdf_3D(SMAP_CRN.to_numpy(copy= True), 2, 3, 4)/N    ## Joint entropy of H(tbh,tbv,Teff)
    ##4 variate entroy
    H_h_v_temp_insitu = pdf.pdf_4D(SMAP_CRN.to_numpy(copy= True), 0, 2, 3, 4)/N   ## Joint entropy of H(tbh,tbv,
                                                                    ##             in situ,Teff)
    ###MI 2 variables
    MI_dca_insitu = H_DCA + H_insitu - H_dca_insitu     ## I(DCA;in situ) = H(in situ) + H(DCA) - H(DCA,in situ)    
    MI_h_dca = H_Tbh + H_DCA - H_h_dca                  ## I(tbh;DCA) = H(tbh) + H(DCA) - H(tbh,DCA)
    MI_v_dca = H_Tbv + H_DCA - H_v_dca                  ## I(tbv;DCA) = H(tbv) + H(DCA) - H(tbv,DCA)
    MI_h_v = H_Tbh + H_Tbv - H_h_v                      ## I(tbh;tbv) = H(tbv) + H(tbh) - H(tbv,tbh)
    ###MI 3 variables 
    MI_hv_dca = H_DCA + H_h_v - H_h_v_dca               ## I(tbh,tbv;DCA) = H(DCA) + H(tbh,tbv) - H(tbh,tbv,DCA)
    ###MI 4 variables 
    MI_hv_temp_insitu = H_insitu + H_h_v_temp - H_h_v_temp_insitu   ## I(tbh,tbv,Teff;in situ) = H(in situ) + H(tbh,tbv,Teff) - H(tbh,tbv,Teff,in situ)    
    ######partial information decompostion 
    ## Here Xtar --> DCA; 
    ##      Xs1  --> TBh
    ##      Xs2  --> TBv
    ##II = I(Xs1;Xtar|Xs2) - I(Xs1;Xtar)
    ##II = S - R
    ##I(tbh,tbv;DCA) = Uh + Uv + S + R
    ##I(tbh;DCA) = Uh + R
    ##I(tbv;DCA) = Uv + R
    ##II = S - R = I(tbh,tbv;DCA) - I(tbh;DCA) - I(tbv;DCA)
    II = MI_hv_dca - MI_h_dca - MI_v_dca       ##The interaction information 
    II2 = H_h_dca + H_v_dca + H_h_v - H_Tbh - H_Tbv - H_DCA - H_h_v_dca
    assert round(II,4) == round(II2,4), "we are in trouble!!!"
    
    ##RMMI = min[I(Xs1;Xtar), I(Xs2;Xtar)]
    ##     = min[I(XTBh;Xtar), I(XTBv;Xtar)]
    R_MMI = min(MI_h_dca,MI_v_dca)

    ##Rmin = max(0, -II)
    R_min = max(0,-II)
    
    ##Source dependcy Is = I(Xs1;Xs2)/min[H(Xs1), H(Xs2))]
    ##                   = I(TBh;TBv)/min[H(TBh), H(TBv)]
    I_s = MI_h_v/min(H_Tbh, H_Tbv)
    ##Rescaled redundancy Rs = Rmin + Is(RMMI - Rmin)
    R_s = R_min + I_s*(R_MMI - R_min)
    
    ## Uh = I(TBh;DCA) - R_s
    U_h = MI_h_dca - R_s
    
    ## Uv = I(TBv;DCA) - R_s
    U_v = MI_v_dca - R_s
    
    ## S = II + R
    S = II + R_s

    metrics.append([i, siteLat, siteLon, sLC,
                    H_DCA, H_insitu, H_Tbh, H_Tbv, H_Teff,
                    MI_dca_insitu, MI_hv_dca, MI_hv_temp_insitu, 
                    U_h, U_v, S, R_s, corrSoilMoisture, corrHandV, 
                    mean_LAI, std_LAI, numPoints, rangeHV])
        
DATA = pd.DataFrame(metrics, columns = ['siteName','siteLat', 'siteLon', 'siteLandCover',
                                        'H_DCA', 'H_insitu', 'H_Tbh', 'H_Tbv', 'H_Teff', 
                                        'MI_dca_insitu','MI_hv_dca', 'MI_hv_temp_insitu',
                                        'U_h', 'U_v', 'S', 'R_s', 'corrSM', 'corrHV', 
                                        'mean_LAI', 'std_LAI', 'numPoints', 'range'])

DATA.to_csv("C:/Users/libon/Box/USCRM_SMAP_sm/HESS_06_08_Loop_SMAP_Final_check_dropdup.csv")
              
'''This script has been double checked on 2021 06 08 14:41:00 PST'''
'''This script has been double checked on 2021 06 25 00:28:00 PST'''


