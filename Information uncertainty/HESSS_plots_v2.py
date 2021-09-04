# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 11:28:40 2020

@author: libon
"""
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt 
from scipy.stats import pearsonr
import seaborn as sns
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from texttable import Texttable
plt.style.use('default')
######folder locations 

Savefig_path = 'C:/Users/libon/Box/USCRM_SMAP_sm/HESS_figures/'
Box = 'C:/Users/libon/Box/' 
UScrn_SMAP = Box + 'USCRM_SMAP_sm/HESS_06_08_Loop_SMAP_Final_check_dropdup.csv'
plt.rcParams["figure.dpi"] = 300




pid_df = pd.read_csv(UScrn_SMAP, index_col=0)
pid_df['long_rank'] = pid_df['siteLon'].rank(ascending=True)

color_df = {'Grasslands':'g','Shrublands':'blue',
            'Croplands':'r','Mixed':'black', 
            'Savannas':'cyan'}
custom_dict = {'Shrublands': 0, 'Grasslands': 1, 
               'Croplands': 2, 'Savannas':3, 
               'Mixed':4,}
spatialMap = pid_df.copy() 
spatialMap['color'] = spatialMap['siteLandCover'].map(color_df) 
spatialMap = spatialMap.sort_values(by=['siteLandCover'], 
                                    key=lambda x: x.map(custom_dict),
                                    ignore_index=True)

def HESSFigure1(data):
    data = data.copy()
    df_loc = data[['siteLon','siteLat', 'siteLandCover', 'color']]
    
    canada_east = -80
    canada_west = -120
    canada_north = 49
    canada_south = 23.5
    standard_parallels = (49, 77)
    central_longitude = -104
    fig1,ax1 = plt.subplots()
    ax1 = plt.axes(projection=ccrs.LambertConformal(central_longitude=central_longitude,
                                                      standard_parallels=standard_parallels))
    
    ax1.set_extent([canada_west, canada_east, canada_south, canada_north],ccrs.Geodetic())
    # ax.stock_img()
    # ax.gridlines(linestyle='--')
    for i in df_loc.to_numpy(copy=True):
        ax1.scatter(i[0],i[1], transform=ccrs.PlateCarree(),
                s = 22, edgecolors = i[-1], facecolor='none', label = i[-2])#
    # Step-by-step
    ax11 = plt.gca()                      
    handle_labels = ax11.get_legend_handles_labels() 

    
    unique_labels = {l:h for h,l in zip(*handle_labels)} 
    
    label_list = [*zip(*unique_labels.items())]
    # print(unique_labels)
    get_list = label_list[::-1]                        
    ax1.legend(*get_list, fontsize = 7, loc = 4)
    ax1.add_feature(cfeature.BORDERS,linewidth = 0.5)
    ax1.add_feature(cfeature.COASTLINE,linewidth = 0.5)
    ax1.add_feature(cfeature.STATES, linewidth = 0.3, linestyle = ':')
    # fig1.subplots_adjust(left=0.01, right=0.99, top=1, bottom=0)

                        
    fig1.savefig(Savefig_path + 'HESS_fg1.png',bbox_inches='tight',pad_inches = 0.1)
##########plotting selected VOD and SM at a given loacation 
HESSFigure1(spatialMap)


def getRename(data):
    data = data.copy()
    loc = []
    for i in data['siteName'].to_numpy():
        newname = i.split('_')
        newname.pop(-2)
        newname.pop(0)
        Name = ' '.join(newname) if newname[0] in ['Newton', 'Goodwell', 'Wolf', 'Stillwater']  else ' '.join(newname[:-1])
        loc.append(Name)
    data['newname'] = loc
    return data




def HESSFigure2(data):
    data = data.copy()
    plt.style.use('seaborn')
    data = getRename(data)   
    data = data.sort_values(by = 'siteLon') ##sort by long
    fig2, ax2 = plt.subplots(figsize=(12,7))

    ax2.plot(data['newname'], data['H_Tbh'], label = '$H_{CN}(T_{Bh})$', 
            linestyle='-', marker='s', c= 'blue')
    ax2.plot(data['newname'], data['H_Tbv'], label = '$H_{CN}(T_{Bv})$', 
            linestyle='-', marker='s', c= 'cyan')
    ax2.plot(data['newname'], data['H_insitu'], label = '$H_{CN}(in\ situ)$', 
            linestyle='--', c= 'red', marker = 'D') 
    ax2.plot(data['newname'], data['H_DCA'], label = '$H_{CN}(DCA)$', 
            linestyle='-', marker='v', c= 'green')
    ax2.plot(data['newname'], data['H_Teff'], label = '$H_{CN}(T_{eff})$', 
            linestyle='-', marker='o', c= 'black')
    ax2.set_ylim([0.2, 0.45])
    ax2.set_xlabel('Study sites', fontsize = 16)
    ax2.tick_params(axis='both', labelsize = 13)
    ax2.set_ylabel('Entropy',fontsize = 16)
    ax2.tick_params(axis='x', labelrotation = 90)
    ax2.legend(fontsize = 13, loc = 3, columnspacing = 0.15)
    fig2.savefig(Savefig_path + 'HESS_fg2v2.png',bbox_inches='tight',pad_inches = 0.1) 
        
HESSFigure2(spatialMap)


def HESSFigure3(data):
    data =data.copy()
    data = getRename(data) 
    data = data.sort_values(by = 'siteLon')
    fig3, ax3 = plt.subplots(figsize=(12,7))
    ax3.plot(data['newname'], data['MI_hv_temp_insitu'], label = '$I(T_{Bh},T_{Bv},T_{eff};in\ situ)$',
            linestyle='-', marker='D', c = 'g')
    ax3.plot(data['newname'], data['MI_hv_dca'], label = '$I(T_{Bh},T_{Bv};DCA)$',
            linestyle='-', marker='o', c = 'b')
    ax3.plot(data['newname'], data['MI_dca_insitu'], label = '$I(DCA;in\ situ)$',
            linestyle='-', marker='s', c = 'r')

    # ax.set_ylim([0, 0.67])
    ax3.set_xlabel('Study sites', fontsize = 16)
    ax3.tick_params(axis='both', labelsize = 13)
    ax3.set_ylabel('Mutual information',fontsize = 16)
    ax3.tick_params(axis='x', labelrotation = 90)
    ax3.legend(fontsize = 13, loc = 2)
    ax3.set_ylim([0, 0.28])
    fig3.savefig(Savefig_path + 'HESS_fg3v2.png',
                  bbox_inches='tight',pad_inches = 0.1)  
        
HESSFigure3(spatialMap)





def HESSFigure4(data):
    data = data.copy()
    fig4, ax4 = plt.subplots(1,2, figsize=(16,8))
    data = data[['H_Tbh','H_Tbv','H_insitu', 'H_DCA', 'H_Teff',
                 'MI_hv_temp_insitu', 'MI_hv_dca', 'MI_dca_insitu']]   
                      
    print('The following is the HESS Figure4a!')
    tableFigure4a = Texttable()
    R_Tbh,  pR_h = pearsonr(data['H_insitu'], data['H_Tbh'])
    R_Tbv,  pR_v = pearsonr(data['H_insitu'], data['H_Tbv'])
    R_Teff, pR_Teff = pearsonr(data['H_insitu'], data['H_Teff']) 
    R_Dca,  pR_Dca = pearsonr(data['H_insitu'], data['H_DCA'])

    SpR_h = '<' if pR_h < 0.05 else '>'
    SpR_v = '<' if pR_v < 0.05 else '>'
    SpR_Teff = '<' if pR_Teff < 0.05 else '>'
    SpR_Dca = '<' if pR_Dca < 0.05 else '>'
    
    tableFigure4a.add_rows([['R(insitu,Tbh), pvalue', 'R(insitu,Tbv), pvalue',
                      'R(insitu,Teff), pvalue', 'R(insitu,DCA), pvalue'],
                    [(R_Tbh, pR_h), (R_Tbv, pR_v), 
                      (R_Teff, pR_Teff), (R_Dca,pR_Dca)]])
    print(tableFigure4a.draw())
    print('\n')
    print('The following is HESS Figure4b')
    sns.regplot(x = data['H_insitu'], y = data['H_Tbh'], ax=ax4[0], 
            line_kws={"linestyle": '-.', 'color':'blue'}, ci = None, 
            scatter_kws={'s': 70,  'color':'blue'}, marker = 's', 
            label = '$H_{CN}(T_{Bh})' + '\ \ r = {},\ '.format(round(R_Tbh,2)) + 'p' + 
            SpR_h + ' 0.05$')
    
    sns.regplot(x = data['H_insitu'], y = data['H_Tbv'], ax=ax4[0], 
            line_kws={"linestyle": '-.', 'color':'cyan'}, ci = None, 
            scatter_kws={'s': 70,  'color':'cyan'}, marker = 'D', 
            label = '$H_{CN}(T_{Bv})' + '\ \ r = {},\ '.format(round(R_Tbv,2)) + 'p' + 
            SpR_v + ' 0.05$')
    sns.regplot(x =data['H_insitu'], y = data['H_DCA'], ax=ax4[0], 
            line_kws={"linestyle": '-.', 'color':'green'}, ci = None, 
            scatter_kws={'s': 70,  'color':'green'}, marker = 'v', 
            label = '$H_{CN}(DCA)' + '\ \ r = {},\ '.format(round(R_Dca,2)) + 'p' + 
            SpR_Dca + ' 0.05$')
    sns.regplot(x = data['H_insitu'], y = data['H_Teff'], ax=ax4[0], 
            line_kws={"linestyle": '-.', 'color':'black'}, ci = None, 
            scatter_kws={'s': 70,  'color':'black'}, marker = 'o', 
            label = '$H_{CN}(T_{eff})' + '\ \ r = {},\ '.format(round(R_Teff,2)) + 'p' + 
            SpR_Teff + ' 0.05$')
    
    ax4[0].set_xlim(0.26, 0.411)
    ax4[0].set_ylim(0.275, 0.46)

    ax4[1].set_xlim(0.26, 0.411)
    ax4[1].set_ylim(0, 0.27)
    ####################
    tableFigure4b = Texttable()
    r_itt, pitt = pearsonr(data['H_insitu'], data['MI_hv_temp_insitu'])
    r_idcahv, pidcahv = pearsonr(data['H_insitu'], data['MI_hv_dca'])
    r_idcainsitu, pidcainsitu = pearsonr(data['H_insitu'], data['MI_dca_insitu']) ###not significant
    
    Sig_pitt = '<' if pitt < 0.05 else '>'
    Sig_pidcahv = '<' if pidcahv < 0.05 else '>'
    Sig_pidcainsitu = '<' if pidcainsitu < 0.05 else '>'
    tableFigure4b.add_rows([['R(insitu,I(H,V,Teff;insitu)), pvalue',
                              'R(insitu,I(DCA;H,V)), pvalue',
                              'R(insitu,I(DCA;insitu)), pvalue'],
                    [(r_itt, pitt), (r_idcahv,pidcahv), (r_idcainsitu,pidcainsitu)]])
    print(tableFigure4b.draw())
    print('\n')
    
    sns.regplot(x = data['H_insitu'], y = data['MI_hv_temp_insitu'], ax=ax4[1], 
            line_kws={"linestyle": '-.', 'color':'g'}, ci = None, 
            scatter_kws={'s': 70,  'color':'g'}, marker = 'D', 
            label = '$I(T_{Bh},T_{Bv},T_{eff};in\ situ)' + '\ \ r = {},\ '.format(round(r_itt,2)) + 'p' + 
            Sig_pitt + ' 0.05$')
    
    sns.regplot(x = data['H_insitu'], y = data['MI_hv_dca'], ax=ax4[1], 
            line_kws={"linestyle": '-.', 'color':'b'}, ci = None, 
            scatter_kws={'s': 70,  'color':'b'}, marker = 'o', 
            label = '$I(T_{Bh},T_{Bv};DCA)' + '\ \ r = {:.2f},\ '.format(round(r_idcahv,2)) + 'p' + 
            Sig_pidcahv + ' 0.05$')
    
    sns.regplot(x= data['H_insitu'], y = data['MI_dca_insitu'], ax=ax4[1], 
            line_kws={"linestyle": '-.', 'color':'r'}, ci = None, 
            scatter_kws={'s': 70,  'color':'r'}, marker = 's', 
            label = '$I(DCA;in\ situ)' + '\ \ r = {:.2f},\ '.format(round(r_idcainsitu,2)) + 'p' + 
            Sig_pidcainsitu + ' 0.05$')
               
    ax4[0].tick_params(axis='both', labelsize = 15)
    ax4[1].tick_params(axis='both', labelsize = 15)
    ax4[0].legend(fontsize = 13, loc = 2)          
    ax4[1].legend(fontsize = 13, loc = 2)    
    ax4[0].set_xlabel('Entropy of $in\ situ$, $H_{CN}(in\ situ)$', fontsize = 17)
    ax4[0].set_ylabel('Entropy',fontsize = 17)
    ax4[1].set_xlabel('Entropy of $in\ situ$, $H_{CN}(in\ situ)$', fontsize = 17)    
    ax4[1].set_ylabel('Mutual information',fontsize = 17)
 
 

    ax4[0].text(0.85,0.9,'a', fontsize = 18, transform = ax4[0].transAxes)
    ax4[1].text(0.85,0.9,'b', fontsize = 18, transform = ax4[1].transAxes)
    plt.subplots_adjust(wspace=0.2)
    fig4.savefig(Savefig_path + 'HESS_fg4v2.png', bbox_inches='tight',pad_inches = 0.1) 
HESSFigure4(spatialMap)




######Scatter plot of Model uncertainty vs Pearson r
def HESSFigure5(data):
    data = data.copy()
    fig5, ax5 = plt.subplots(3, 1, figsize=(11,24))
    data['total_error'] = data['H_insitu']- data['MI_dca_insitu']
    data['model_error'] = data['MI_hv_temp_insitu'] - data['MI_dca_insitu']
    data['random_error'] = data['total_error'] - data['model_error']

    data = data[['total_error', 'model_error', 'random_error','corrSM', 
                  'color','siteLandCover']]   
    
    min_total, max_total = min(data['total_error']), max(data['total_error'])
    min_model, max_model =  min(data['model_error']), max(data['model_error'])
    min_random, max_random =  min(data['random_error']), max(data['random_error'])
    
    R_mod = pearsonr(data['model_error'],  data['corrSM'])
    R_rnd = pearsonr(data['random_error'], data['corrSM'])
    R_tot = pearsonr(data['total_error'],  data['corrSM'])
    
    corrs = [R_tot,R_mod, R_rnd]
    
    tableFigure5a = Texttable()
    print("The following print informational uncertainty (Figure5) vs. pearsonr(DCA;insitu)")
    tableFigure5a.add_rows([['R(TotalError,R(insitu;DCA)), pvalue',
                             'R(ModelError,R(insitu;DCA)), pvalue',
                             'R(RadomError,R(insitu;DCA)), pvalue'],
                            [R_tot, R_mod, R_rnd]])
    print(tableFigure5a.draw())
    print('\n')
    
    sns.regplot(x = data['total_error'], y = data['corrSM'],ax=ax5[0], 
            line_kws={"linestyle": '-.', 'color':'black'},
            scatter=False)
    sns.regplot(x = data['model_error'], y = data['corrSM'],ax=ax5[1], 
            line_kws={"linestyle": '-.', 'color':'black'},
            scatter=False)
    sns.regplot(x= data['random_error'], y = data['corrSM'],ax=ax5[2], 
            line_kws={"linestyle": '-.', 'color':'black'},
            scatter=False)
    scatterSize = 60
    for i in data.to_numpy():
        ax5[0].scatter(i[0], i[3], c = i[-2], label = i[-1], s= scatterSize)
        ax5[1].scatter(i[1], i[3], c = i[-2], label = i[-1], s= scatterSize)
        ax5[2].scatter(i[2], i[3], c = i[-2], label = i[-1], s= scatterSize)
    def getHandleLable(hl):
        hdls  = hl[0]
        lbs = hl[1]
        hOut = []
        lOut = []
        for h, l in zip(hdls, lbs):
            if l not in lOut:
                lOut.append(l)
                hOut.append(h)
        return hOut, lOut
    ax20Handle,ax20Labels = getHandleLable(ax5[0].get_legend_handles_labels())
    ax21Handle,ax21Labels = getHandleLable(ax5[1].get_legend_handles_labels())
    ax22Handle,ax22Labels = getHandleLable(ax5[2].get_legend_handles_labels())
    ax5[0].legend(ax20Handle,ax20Labels,
                  fontsize = 20, loc = 3)  
    ax5[1].legend(ax21Handle,ax21Labels,
                  fontsize = 20, loc = 3),
    ax5[2].legend(ax22Handle,ax22Labels,
                  fontsize = 20, loc = 3)
    [ax5[i].set_ylim(0, 1) for i in range(3)]
    ax5[0].set_xlim(0.995*min_total, max_total*1.005)
    ax5[1].set_xlim(0.96*min_model, max_model*1.007)   
    ax5[2].set_xlim(0.993*min_random, max_random*1.005)
    ax5[0].set_xlabel('SMAP informational total uncertainty, $I_{Tot}$',fontsize = 21)
    ax5[1].set_xlabel('SMAP informational model uncertainty, $I_{Mod}$',fontsize = 21)
    ax5[2].set_xlabel('SMAP informational random uncertainty, $I_{Rnd}$',fontsize = 21)
    [ax5[i].set_ylabel('Pearson correlation between \n SMAP  DCA and $in\ situ$ (unitless)',fontsize = 21) 
      for i in range(3)]
    [ax5[i].tick_params(axis = 'both', labelsize = 18) for i in range(3)]
   
    [ax5[i].text(0.65,0.89, '$ r = {},\ p < {}$'.format(round(corrs[i][0],2),0.05), 
              fontsize = 20,transform = ax5[i].transAxes) for i in np.arange(len(corrs))]
    numberLabels = ['a','b','c']
    [ax5[i].text(0.85,0.75, numberLabels[i], fontsize = 25, transform = ax5[i].transAxes) 
     for i in np.arange(len(numberLabels))]

    plt.subplots_adjust(hspace = 0.16)
    fig5.savefig(Savefig_path + 'HESS_fg5v2.png',
                  bbox_inches='tight',pad_inches = 0.1) 
HESSFigure5(spatialMap)



#########PID Bar plot
def HESSFigure6(data):
    # fig4, ax4 = plt.subplots()
    data = data.copy()
    data = getRename(data)
    data = data.sort_values(by = 'siteLon')
    colorDict = {'R_s':'r', 'S':'g', 'U_h' :'b', 'U_v':'y'}
    ax6 = data[['R_s','S','U_h','U_v', 'newname']].plot.bar(x = 'newname', 
                                                            y =['R_s','S','U_h','U_v'], 
                                                            figsize=(50,30),fontsize = 50,
                                                            rot = 90, color = colorDict,
                                                            width=0.6) 
    ax6.legend(['Redundant information of $T_{Bh}$ and $T_{Bv}$','Synergistic information of $T_{Bh}$ and $T_{Bv}$',
                'Unique information of $T_{Bh}$','Unique information of $T_{Bv}$'],fontsize=50, loc = 1)
    ax6.set_ylabel('Partial information decomposition components',fontsize=60)
    ax6.set_xlabel('Study sites',fontsize=60)
    plt.savefig(Savefig_path + 'HESS_fg6v2.png', 
                bbox_inches='tight',pad_inches = 0.1) 
HESSFigure6(spatialMap)



####PID plots VS pearson r
def HESSFigure7(data):
    data = data.copy()
    fig7, ax7 = plt.subplots(2,2, figsize = (18,18))  
   
    newlabels = ['Unique information of $T_{Bh}$','Unique information of $T_{Bv}$',
                  'Synergistic information of $T_{Bh}$ and $T_{Bv}$',
                  'Redundant information of $T_{Bh}$ and $T_{Bv}$']
    old_labels = ['U_h','U_v','S','R_s'] 
    renamedData = data.rename(columns= dict(zip(old_labels, newlabels)))
    subplot_num = ['a','b','c','d']
    for i in np.arange(len(newlabels)):
        infoComponent = newlabels[i]
        N_components = renamedData[[infoComponent,'corrSM','color', 'siteLandCover']]
        corr = pearsonr(N_components[infoComponent], N_components['corrSM'])
        if corr[1] < 0.05: 
            ax7[i//2][i%2].text(0.40,0.91,'$ r = {:.2f},\ p < {}$'.format(corr[0],0.05),transform=ax7[i//2][i%2].transAxes,
                                fontsize = 19)
            Line = True 
        else:
            Line = False
        for j in N_components.values:
            ####now plot the information components vs. pearson correlations
            ax7[i//2][i%2].scatter(j[0], j[1], s = 90,
                                      c = j[2], label = j[3])
        handle,label = ax7[i//2][i%2].get_legend_handles_labels() ##handle_label = [(h1...),(l2...)]
        handleOut = []
        labelOut = []
        labelHandle = Texttable()
        labelHandle.add_row(['lable'])
        for h, l in zip(handle, label):
            if l not in labelOut:
                labelHandle.add_row([l])
                labelOut.append(l)
                handleOut.append(h)
        print(labelHandle.draw())
        print()
        ax7[i//2][i%2].legend(handleOut,labelOut,
                  fontsize = 17, loc = 3,
                  bbox_transform = ax7[i//2][i%2].transAxes)
                    
        sns.regplot(x = N_components[infoComponent],
                    y = N_components['corrSM'],ax=ax7[i//2,i%2],
                    line_kws={"linestyle": '-.', 'color':'black'},
                    scatter = False, fit_reg= Line)
   
        ax7[i//2][i%2].text(0.9,0.8,subplot_num[i],transform=ax7[i//2][i%2].transAxes,
                                fontsize = 23) 

        ax7[i//2][i%2].set_ylim(0, 1)
                      
        if i in [1,2]:
            ax7[i//2][i%2].set_xlim(0.85*min(N_components[infoComponent]),1.01*max(N_components[infoComponent]))
        elif i == 0:
            ax7[i//2][i%2].set_xlim(0.6*min(N_components[infoComponent]),1.01*max(N_components[infoComponent]))
        else:
            ax7[i//2][i%2].set_xlim(0.94*min(N_components[infoComponent]),1.01*max(N_components[infoComponent]))
        if i in [0,2]:
            ax7[i//2][i%2].set_ylabel('Pearson correlation between \n SMAP  DCA and $in\ situ$ (unitless)'
                                      ,fontsize = 22)
        else:
            ax7[i//2][i%2].yaxis.label.set_visible(False)
        ax7[i//2][i%2].set_xlabel(infoComponent,fontsize = 22)            
        ax7[i//2][i%2].tick_params(axis = 'both', labelsize = 18)


    plt.subplots_adjust(wspace=0.13, hspace=0.15)
    fig7.savefig(Savefig_path + 'HESS_fg7v2.png',
                  bbox_inches='tight',pad_inches = 0.1) 
HESSFigure7(spatialMap)




'''This script was double checked on 2021 06 08 16:30:00'''


