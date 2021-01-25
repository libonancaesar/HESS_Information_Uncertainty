# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 11:28:40 2020

@author: libon
"""
import os 
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt 
from scipy.stats import pearsonr, ttest_ind, ttest_rel
from scipy import stats
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
plt.style.use('default')
######folder locations 
Box = 'C:/Users/libon/Box/' 
UScrn_SMAP = Box + 'USCRM_SMAP_sm/HESS_PID_MM.csv'
plt.rcParams["figure.dpi"] = 300
Savefig_path = 'C:/Users/libon/Dropbox/Bonan_data and scripts/HESS_2020/HESS_2020_figure/'
######
pid_df = pd.read_csv(UScrn_SMAP, index_col=0)
pid_df['long_rank'] = pid_df['long'].rank(ascending=True)
color_df = {'Grasslands':'g','Shrublands':'blue','Croplands':'r','Mixed':'black'}
custom_dict = {'Shrublands': 0, 'Grasslands': 1, 'Croplands': 2, 'Mixed': 3} 
pid_df['color'] = pid_df['landcover'].map(color_df) 
pid_df = pid_df.sort_values(by=['landcover'], key=lambda x: x.map(custom_dict), ignore_index=True)

def getFigure1(data):
    df_loc = data[['long','lat', 'landcover', 'color']]
    
    canada_east = -80
    canada_west = -120
    canada_north = 49
    canada_south = 23.5
    standard_parallels = (49, 77)
    central_longitude = -104
    fig,ax = plt.subplots()
    ax = plt.axes(projection=ccrs.LambertConformal(central_longitude=central_longitude,
                                                      standard_parallels=standard_parallels))
    
    ax.set_extent([canada_west, canada_east, canada_south, canada_north],ccrs.Geodetic())
    # ax.stock_img()
    # ax.gridlines(linestyle='--')
    for i in df_loc.values:
        ax.scatter(i[0],i[1], transform=ccrs.PlateCarree(),
                s = 22, edgecolors = i[-1], facecolor='none', label = i[-2])#
    # Step-by-step
    ax1 = plt.gca()                      
    handle_labels = ax1.get_legend_handles_labels() 
    unique_labels = {l:h for h,l in zip(*handle_labels)}                 
    label_list = [*zip(*unique_labels.items())]             
    get_list = label_list[::-1]                        
    ax.legend(*get_list, fontsize = 7, loc = 4)
    ax.add_feature(cfeature.BORDERS,linewidth = 0.5)
    ax.add_feature(cfeature.COASTLINE,linewidth = 0.5)
    ax.add_feature(cfeature.STATES, linewidth = 0.3, linestyle = ':')
    fig.savefig(Savefig_path + 'uscrn_study_map.png')
##########plotting selected VOD and SM at a given loacation 
getFigure1(pid_df)
plt.style.use('seaborn')


def getRename(data):
    loc = []
    for i in data['sitename'].values:
        newname = i.split('_')
        newname.pop(-2)
        newname.pop(0)
        Name = ' '.join(newname) if newname[0] in ['Newton', 'Goodwell', 'Wolf']  else ' '.join(newname[:-1])
        loc.append(Name)
    data['newname'] = loc
    return data




def getFigure3(data):
    data = getRename(data)    
    fig, ax = plt.subplots(figsize=(12,7))

    ax.plot(data['newname'], data['H_h_36'], label = '$H_{CN}(T_{Bh})$', 
            linestyle='-', marker='s', c= 'blue')
    ax.plot(data['newname'], data['H_v_36'], label = '$H_{CN}(T_{Bv})$', 
            linestyle='-', marker='s', c= 'cyan')
    ax.plot(data['newname'], data['H_in_situ'], label = '$H_{CN}(in\ situ)$', 
            linestyle='--', c= 'red', marker = 'D') 
    ax.plot(data['newname'], data['H_dca36'], label = '$H_{CN}(DCA)$', 
            linestyle='-', marker='v', c= 'green')
    ax.plot(data['newname'], data['H_temp36'], label = '$H_{CN}(T_{eff})$', 
            linestyle='-', marker='o', c= 'black')
    ax.set_ylim([0.2, 0.45])
    ax.set_xlabel('Study sites', fontsize = 16)
    ax.tick_params(axis='both', labelsize = 13)
    ax.set_ylabel('Entropy',fontsize = 16)
    ax.tick_params(axis='x', labelrotation = 90)
    ax.legend(fontsize = 13, loc = 3, columnspacing = 0.15)
    data = data.reindex(np.arange(len(data)))  
    
    # for index in data.index:
    #     ax.get_xticklabels()[index].set_color(data['color'][index])
        
        
getFigure3(pid_df)


def getFigure4(data):
    data = getRename(data)    
    fig, ax = plt.subplots(figsize=(12,7))
    ax.plot(data['newname'], data['I_hv_temp_crn36'], label = '$I(T_{Bh},T_{Bv},T_{eff};in\ situ)$',
            linestyle='-', marker='D', c = 'g')
    ax.plot(data['newname'], data['I_hv_dca_36'], label = '$I(T_{Bh},T_{Bv};DCA)$',
            linestyle='-', marker='o', c = 'b')
    ax.plot(data['newname'], data['I_crn_dca36'], label = '$I(DCA;in\ situ)$',
            linestyle='-', marker='s', c = 'r')

    # ax.set_ylim([0, 0.67])
    ax.set_xlabel('Study sites', fontsize = 16)
    ax.tick_params(axis='both', labelsize = 13)
    ax.set_ylabel('Mutual information',fontsize = 16)
    ax.tick_params(axis='x', labelrotation = 90)
    ax.legend(fontsize = 13, loc = 2)
    ax.set_ylim([0, 0.28])
    data = data.reindex(np.arange(len(data)))  
    
    # for index in data.index:
    #     ax.get_xticklabels()[index].set_color(data['color'][index])
        
        
getFigure4(pid_df)





def getFigure5(data):
    fig1, ax1 = plt.subplots(1,2, figsize=(16,8))
    used_Data = data[['H_h_36','H_v_36','H_in_situ', 'H_dca36', 'H_temp36',
                      'I_hv_temp_crn36', 'I_hv_dca_36', 'I_crn_dca36', 'landcover', 'color']]    
    
    
    hr = pearsonr(data['H_in_situ'], data['H_h_36'])
    vr = pearsonr(data['H_in_situ'], data['H_v_36'])
    tr = pearsonr(data['H_in_situ'], data['H_temp36']) ###not significant
    dcar = pearsonr(data['H_in_situ'], data['H_dca36'])
    print(hr, vr, tr, dcar)
    
    sns.regplot(data['H_in_situ'], data['H_h_36'], ax=ax1[0], 
            line_kws={"linestyle": '-.', 'color':'blue'}, ci = None, 
            scatter_kws={'s': 70,  'color':'blue'}, marker = 's', 
            label = '$H_{CN}(T_{Bh})' + '\ \ r = {:.2f},\ p < {}$'.format(hr[0], 0.05))
    
    sns.regplot(data['H_in_situ'], data['H_v_36'], ax=ax1[0], 
            line_kws={"linestyle": '-.', 'color':'cyan'}, ci = None, 
            scatter_kws={'s': 70,  'color':'cyan'}, marker = 'D', 
            label = '$H_{CN}(T_{Bv})' + '\ \ r = {:.2f},\ p < {}$'.format(vr[0], 0.05))
    
    sns.regplot(data['H_in_situ'], data['H_dca36'], ax=ax1[0], 
            line_kws={"linestyle": '-.', 'color':'green'}, ci = None, 
            scatter_kws={'s': 70,  'color':'green'}, marker = 'v', 
            label = '$H_{CN}(DCA)' + '\ \ r = {:.2f},\ p < {}$'.format(dcar[0], 0.05))
    
    sns.regplot(data['H_in_situ'], data['H_temp36'], ax=ax1[0], 
            line_kws={"linestyle": '-.', 'color':'black'}, ci = None, 
            scatter_kws={'s': 70,  'color':'black'}, marker = 'o', 
            label = '$H_{CN}(T_{eff})' + '\ \ r = {:.2f},\ p > {}$'.format(tr[0], 0.05))
    
    
    ax1[0].set_xlim(0.26, 0.40)
    ax1[0].set_ylim(0.275, 0.46)

    ax1[1].set_xlim(0.26, 0.40)
    ax1[1].set_ylim(0, 0.25)
    
    itttr = pearsonr(data['H_in_situ'], data['I_hv_temp_crn36'])
    ittr = pearsonr(data['H_in_situ'], data['I_hv_dca_36'])
    idcar = pearsonr(data['H_in_situ'], data['I_crn_dca36']) ###not significant
    
    print(itttr, ittr, idcar)
    sns.regplot(data['H_in_situ'], data['I_hv_temp_crn36'], ax=ax1[1], 
            line_kws={"linestyle": '-.', 'color':'g'}, ci = None, 
            scatter_kws={'s': 70,  'color':'g'}, marker = 'D', 
            label = '$I(T_{Bh},T_{Bv},T_{eff};in\ situ)' + '\ \ r = {:.2f},\ p < {}$'.format(itttr[0], 0.05))
    
    sns.regplot(data['H_in_situ'], data['I_hv_dca_36'], ax=ax1[1], 
            line_kws={"linestyle": '-.', 'color':'b'}, ci = None, 
            scatter_kws={'s': 70,  'color':'b'}, marker = 'o', 
            label = '$I(T_{Bh},T_{Bv};DCA)' + '\ \ r = {:.2f},\ p > {}$'.format(ittr[0], 0.05))
               
    sns.regplot(data['H_in_situ'], data['I_crn_dca36'], ax=ax1[1], 
            line_kws={"linestyle": '-.', 'color':'r'}, ci = None, 
            scatter_kws={'s': 70,  'color':'r'}, marker = 's', 
            label = '$I(DCA;In\ situ)' + '\ \ r = {:.2f},\ p > {}$'.format(idcar[0], 0.05))
               
    ax1[0].tick_params(axis='both', labelsize = 15)
    ax1[1].tick_params(axis='both', labelsize = 15)
    ax1[0].legend(fontsize = 13, loc = 2)          
    ax1[1].legend(fontsize = 13, loc = 2)    
    ax1[0].set_xlabel('Entropy of $In\ situ$, $H_{CN}(in\ situ)$', fontsize = 17)
    ax1[0].set_ylabel('Entropy',fontsize = 17)
    ax1[1].set_xlabel('Entropy of $In\ situ$, $H_{CN}(in\ situ)$', fontsize = 17)    
    ax1[1].set_ylabel('Mutual information',fontsize = 17)
 
 

    ax1[0].text(0.9,0.9,'a', fontsize = 18, transform = ax1[0].transAxes)
    ax1[1].text(0.9,0.9,'b', fontsize = 18, transform = ax1[1].transAxes)
    plt.subplots_adjust(wspace=0.2)
getFigure5(pid_df)








#####Uncertainty components 
def getFigure000(data):
    Error = 'I_hv_temp_crn' 
    # dataLabelx = r'$\frac{I(h,v,T;In\ situ)}{H(In\ situ)}$' if var != 'hv' else r'$\frac{I(h,v;In\ situ)}{H(In\ situ)}$'
    dataLabelx =  'Total uncertainty'
    dataLabely = 'Model uncertainty'
    
    
    fig1, ax1 = plt.subplots(figsize=(9,9))
    data['model_error36'] = data[Error + '36'] - data['I_crn_dca36']
    data['model_error9'] = data[Error + '9'] - data['I_crn_dca9']
    
    data['total_error36'] = data['H_in_situ'] - data['I_crn_dca36']
    data['total_error9'] = data['H_in_situ'] - data['I_crn_dca9']
    
    
    data['model_error36_hv'] = data['I_hv_crn36'] - data['I_crn_dca36']
    data['model_error9_hv'] = data['I_hv_crn9'] - data['I_crn_dca9']
    
    lim = max(max(data['total_error36']), max(data['total_error9'])) + 0.2
    t_id = ttest_ind(data.total_error36, data.total_error9)
    t_id_1 = ttest_ind(data.model_error36, data.model_error9)
    print(t_id, t_id_1) ##fail to reject 
    
    
    used_Data = data[['model_error36', 'model_error9', 
                    'total_error36','total_error9', 
                    'model_error36_hv','model_error9_hv',
                    'color', 'landcover']]    

    ax1.scatter(data['model_error36'], data['model_error9'], zorder = 2,
                    s = 90, marker = 'd', label = 'Informational model uncertainty, $I_{Mod}$')
    ax1.scatter(data['total_error36'], data['total_error9'], zorder = 2,
                    s = 90, marker = 'o', label = 'Informational total uncertainty, $I_{Tot}$')
    ax1.set_xlabel('SMAP informational uncertainty at 36km', fontsize=18)
    ax1.set_ylabel('SMAP informational uncertainty at 9km',fontsize=18)   
    ax1.plot([0, lim], [0, lim], c = 'black', linewidth = 2,
                alpha = 0.7, zorder = 1)

    ax1.tick_params(axis = 'both', labelsize = 16)
    # handles0,labels0 = ax1[0].get_legend_handles_labels()
    # handles0 = [handles0[1],handles0[2],handles0[0]]
    # labels0 = [labels0[1], labels0[2], labels0[0]]

    ax1.legend(loc = 2,fontsize = 17, ncol = 1)
    ax1.set_ylim(0, lim)
    ax1.set_xlim(0, lim)
    fig2, ax2 = plt.subplots(figsize=(9,9))

    for i in used_Data.values:
        ax2.scatter(i[2], i[0], c = i[-2], label = i[-1] , s= 90)
        ax2.scatter(i[2], i[4], marker = 'o',
                        edgecolors = i[-2], facecolor = 'none',label = i[-1]+ ' (exclude T)',
                        s= 90, linewidth = 1)
                   
 
    ax2.set_xlabel('SMAP informational total  uncertainty, $I_{Tot}$', fontsize=18)
    ax2.set_ylabel('SMAP informational model  uncertainty, $I_{Mod}$',fontsize=18)
 
    ax2.set_ylim(-0.02, lim)
    ax2.set_xlim(-0.02, lim)
    ax2.tick_params(axis = 'both', labelsize = 16)
    ax2.plot([-0.02,lim],[-0.02,lim], c = 'black', linewidth = 2, 
                alpha = 0.7, zorder = 1)
    handles,labels = ax2.get_legend_handles_labels()
    handles1 = [handles[i] for i in [0,10,-10,-2, 1, 11,-9, -1]]
    labels1 = [labels[i] for i in [0,10,-10,-2, 1, 11,-9, -1]]
    ax2.legend(handles1, labels1, fontsize = 17, 
                  loc = 2, ncol = 2, columnspacing= 0.1) 
    # a = [*zip(*{l:h for h,l in zip(*ax1[1].get_legend_handles_labels())}.items())][::-1]
    # ax1[1].legend(*[*zip(*{l:h for h,l in zip(*ax1[1].get_legend_handles_labels())}.items())][::-1],
    #               fontsize = 17, loc = 2, ncol = 2)    
    # print(a)
    # ax1.text(0.9,0.8,'a', fontsize = 25, transform = ax1.transAxes)
    # ax2.text(0.9,0.8,'b', fontsize = 25, transform = ax2.transAxes)
getFigure000(pid_df)



######Scatter plot of Model uncertainty vs Pearson r
def getFigure4(data):
    fig2, ax2 = plt.subplots(3, 1, figsize=(11,24))
    data['total_error'] = data['H_in_situ']- data['I_crn_dca36']
    data['model_error'] = data['I_hv_temp_crn36'] - data['I_crn_dca36']
    data['random_error'] = data['total_error'] - data['model_error']
    
    data = data[['total_error', 'model_error', 'random_error','Pear_R36', 
                 'RMSE_36','color','landcover']]   
    min_total, max_total = min(data.total_error), max(data.total_error)
    min_model, max_model =  min(data.model_error), max(data.model_error)
    min_random, max_random =  min(data.random_error), max(data.random_error)
    
    R_mod = pearsonr(data['model_error'],data['Pear_R36'])
    R_rnd = pearsonr(data['random_error'],data['Pear_R36'])
    R_tot = pearsonr(data['total_error'],data['Pear_R36'])
    Rs = [R_tot[0],R_mod[0] ,R_rnd[0]]
    print(data.head())
    sns.regplot(data['total_error'],data['Pear_R36'],ax=ax2[0], 
            line_kws={"linestyle": '-.', 'color':'black'},
            scatter=False)
    sns.regplot(data['model_error'],data['Pear_R36'],ax=ax2[1], 
            line_kws={"linestyle": '-.', 'color':'black'},
            scatter=False)
    sns.regplot(data['random_error'],data['Pear_R36'],ax=ax2[2], 
            line_kws={"linestyle": '-.', 'color':'black'},
            scatter=False)
    for i in data.values:
        print(i)
        ax2[0].scatter(i[0], i[3], c = i[-2], label = i[-1], s= 40)
        ax2[1].scatter(i[1], i[3], c = i[-2], label = i[-1], s= 40)
        ax2[2].scatter(i[2], i[3], c = i[-2], label = i[-1], s= 40)
    ax2[0].legend(*[*zip(*{l:h for h,l in zip(*ax2[0].get_legend_handles_labels())}.items())][::-1],
                  fontsize = 20, loc = 3)  
    ax2[1].legend(*[*zip(*{l:h for h,l in zip(*ax2[1].get_legend_handles_labels())}.items())][::-1],
                  fontsize = 20, loc = 3)  
    ax2[2].legend(*[*zip(*{l:h for h,l in zip(*ax2[2].get_legend_handles_labels())}.items())][::-1],
                  fontsize = 20, loc = 3)
    [ax2[i].set_ylim(0, 1) for i in range(3)]
    ax2[0].set_xlim(0.995*min_total, max_total*1.005)
    ax2[1].set_xlim(0.96*min_model, max_model*1.007)   
    ax2[2].set_xlim(0.993*min_random, max_random*1.005)
    ax2[0].set_xlabel('SMAP informational total uncertainty, $I_{Tot}$',fontsize = 21)
    ax2[1].set_xlabel('SMAP informational model uncertainty, $I_{Mod}$',fontsize = 21)
    ax2[2].set_xlabel('SMAP informational random uncertainty, $I_{Rnd}$',fontsize = 21)
    [ax2[i].set_ylabel('Pearson correlation between \n SMAP  DCA and $in\ situ$ (unitless)',fontsize = 21) 
     for i in range(3)]
    [ax2[i].tick_params(axis = 'both', labelsize = 18) for i in range(3)]
   
    [ax2[i].text(0.65,0.89, '$ r = {:.2f},\ p < {}$'.format(Rs[i],0.05), 
              fontsize = 20,transform = ax2[i].transAxes) for i in np.arange(len(Rs))]
    plt.subplots_adjust(hspace = 0.16)
getFigure4(pid_df)



#########PID Bar plot
def getFigure5(data):
    fig4, ax4 = plt.subplots()
    data = getRename(pid_df)
    ax4 = data[['R_36','S_36','U_h_36','U_v_36']].plot.bar(figsize=(50,30),fontsize = 35,
                                             rot = 90, color=['r', 'g', 'b', 'y'],
                                             width=0.7)

    ax4.set_xticks(np.arange(len(data)))
    ax4.set_xticklabels(data['newname'],fontsize = 50)  
    ax4.tick_params(axis = 'both', labelsize = 45)
    for index in data.index:
        ax4.get_xticklabels()[index].set_color(data['color'][index])

    ax4.legend(['Redundant information of $T_{Bh}$ and $T_{Bv}$','Synergistic information of $T_{Bh}$ and $T_{Bv}$',
                'Unique information of $T_{Bh}$','Unique information of $T_{Bv}$'],fontsize=50, loc = 2)
    ax4.set_ylabel('Partial information decomposition components',fontsize=60)
    ax4.set_xlabel('Study sites',fontsize=60)

getFigure5(pid_df)



####PID plots VS pearson r
def getFigure6(data,metric):
    fig5, ax5 = plt.subplots(2,2, figsize = (18,18))   
    newlabels = ['Unique information of $T_{Bh}$','Unique information of $T_{Bv}$',
                 'Synergistic inforamtion of $T_{Bh}$ and $T_{Bv}$',
                 'Redundant information of $T_{Bh}$ and $T_{Bv}$']
    R = '36' 
    old_labels = ['U_h_' + R,'U_v_' + R,'S_' + R,'R_' + R] 
    subplot_num = ['a','b','c','d']
    for i in np.arange(len(old_labels)):
        
        N_components = data[[old_labels[i]] + [metric, 'color', 'landcover']]
        corr = pearsonr(N_components[old_labels[i]], N_components[metric])
        print(corr)
        if corr[1] < 0.05: 
            ax5[i//2][i%2].text(0.42,0.91,'$ r = {:.2f},\ p < {}$'.format(corr[0],0.05),transform=ax5[i//2][i%2].transAxes,
                                fontsize = 19)
            Line = True 
        else:
            Line = False
        for j in N_components.values:
            ax5[i//2][i%2].scatter(j[0], j[1], s = 50,
                                      c = j[2], label = j[3]) 
        a = [*zip(*{l:h for h,l in zip(*ax5[i//2][i%2].get_legend_handles_labels())}.items())][::-1]
        ax5[i//2][i%2].legend(*[*zip(*{l:h for h,l in zip(*ax5[i//2][i%2].get_legend_handles_labels())}.items())][::-1],
                  fontsize = 17, loc = 3, bbox_to_anchor = (0.518, -0.01),
                  bbox_transform = ax5[i//2][i%2].transAxes)
                    
        sns.regplot(N_components[old_labels[i]],N_components[metric],ax=ax5[i//2,i%2],
                    line_kws={"linestyle": '-.', 'color':'black'}, scatter = False,
                      fit_reg= Line)
        
        ax5[i//2][i%2].text(0.9,0.8,subplot_num[i],transform=ax5[i//2][i%2].transAxes,
                                fontsize = 22) 

        ax5[i//2][i%2].set_ylim(0, 1)
        
        if metric[:6] == 'Pear_R':
            ylabel = 'Pearson correlation between \n SMAP  DCA and $in\ situ$ (unitless)'
            
        else:
            ylabel = '$RMSE\ ' + '(DCA,$In\ situ$)'
                  
        if i in [1,2]:
            ax5[i//2][i%2].set_xlim(0.85*min(N_components[old_labels[i]]),1.02*max(N_components[old_labels[i]]))
        elif i == 0:
            ax5[i//2][i%2].set_xlim(0.6*min(N_components[old_labels[i]]),1.02*max(N_components[old_labels[i]]))
        else:
            ax5[i//2][i%2].set_xlim(0.95*min(N_components[old_labels[i]]),1.02*max(N_components[old_labels[i]]))
        if i in [0,2]:
            ax5[i//2][i%2].set_ylabel(ylabel,fontsize = 22)
        else:
           ax5[i//2][i%2].yaxis.label.set_visible(False)
        ax5[i//2][i%2].set_xlabel(newlabels[i],fontsize = 22)            
        ax5[i//2][i%2].tick_params(axis = 'both', labelsize = 18)


    plt.subplots_adjust(wspace=0.13, hspace=0.15)

getFigure6(pid_df, 'Pear_R36')



#########For table 1

pid_df['total_error36'] = pid_df['H_in_situ'] - pid_df['I_crn_dca36']
pid_df['model_error36'] = pid_df['I_hv_temp_crn36'] - pid_df['I_crn_dca36']
pid_df['random_error36'] = pid_df['total_error36'] - pid_df['model_error36']


tableData = pid_df[['H_in_situ','I_crn_dca36', 'I_hv_temp_crn36',
                    'model_error36', 'random_error36','total_error36', 'landcover']]

tableData['tot.H_in_situ'] = tableData['total_error36']/tableData['H_in_situ']
tableData['rnd.tot'] = tableData['random_error36']/tableData['total_error36']
tableData['mod.tot'] = tableData['model_error36']/tableData['total_error36']

overall = tableData.describe()
Grass =  tableData[tableData['landcover'] == 'Grasslands'].describe()
Shrub =  tableData[tableData['landcover'] == 'Shrublands'].describe()
Croplands =  tableData[tableData['landcover'] == 'Croplands'].describe()
Mixed =  tableData[tableData['landcover'] == 'Mixed'].describe()


###############calculate 
pids = pid_df[['I_hv_dca_36','U_h_36', 'U_v_36',
                    'R_36', 'S_36', 'landcover']]

pids['Uh%'] = pids['U_h_36']/pids['I_hv_dca_36']
pids['Uv%'] = pids['U_v_36']/pids['I_hv_dca_36']
pids['R%'] = pids['R_36']/pids['I_hv_dca_36']
pids['S%'] = pids['S_36']/pids['I_hv_dca_36']

overall_pids = pids.describe()
Grass_pids =  pids[pids['landcover'] == 'Grasslands'].describe()
Shrub_pids =  pids[pids['landcover'] == 'Shrublands'].describe()
Croplands_pids =  pids[pids['landcover'] == 'Croplands'].describe()
Mixed_pids =  pids[pids['landcover'] == 'Mixed'].describe()


#########################################################


if False:

    ####PID plots VS pearson r
    def getFigure6666(data,metric):
        fig5, ax5 = plt.subplots(2,2, figsize = (18,18))   
        newlabels = ['Unique information of $T_{Bh}$','Unique information of $T_{Bv}$',
                     'Synergistic inforamtion of $T_{Bh}$ and $T_{Bv}$',
                     'Redundant information of $T_{Bh}$ and $T_{Bv}$']
        R = '36' 
        old_labels = ['U_h_' + R,'U_v_' + R,'S_' + R,'R_' + R] 
        subplot_num = ['a','b','c','d']
        for i in np.arange(len(old_labels)):
            
            N_components = data[[old_labels[i]] + [metric, 'color', 'landcover']]
            corr = pearsonr(N_components[old_labels[i]], N_components[metric])
            print(corr)
            if corr[1] < 0.05: 
                ax5[i//2][i%2].text(0.42,0.91,'$ r = {:.2f},\ p < {}$'.format(corr[0],0.05),transform=ax5[i//2][i%2].transAxes,
                                    fontsize = 19)
                Line = True 
            else:
                Line = False
            for j in N_components.values:
                ax5[i//2][i%2].scatter(j[0], j[1], s = 50,
                                          c = j[2], label = j[3]) 
            a = [*zip(*{l:h for h,l in zip(*ax5[i//2][i%2].get_legend_handles_labels())}.items())][::-1]
            ax5[i//2][i%2].legend(*[*zip(*{l:h for h,l in zip(*ax5[i//2][i%2].get_legend_handles_labels())}.items())][::-1],
                      fontsize = 17, loc = 3, bbox_to_anchor = (0.518, -0.01),
                      bbox_transform = ax5[i//2][i%2].transAxes)
                        
            sns.regplot(N_components[old_labels[i]],N_components[metric],ax=ax5[i//2,i%2],
                        line_kws={"linestyle": '-.', 'color':'black'}, scatter = False,
                          fit_reg= Line)
            
            ax5[i//2][i%2].text(0.9,0.8,subplot_num[i],transform=ax5[i//2][i%2].transAxes,
                                    fontsize = 22) 
    
            # ax5[i//2][i%2].set_ylim(0, 1)
            
            if metric[:6] == 'Pear_R':
                ylabel = 'Pearson correlation between \n SMAP  DCA and $in\ situ$ (unitless)'
                
            else:
                ylabel = 'MAPE'
                      
            # if i in [1,2]:
            #     ax5[i//2][i%2].set_xlim(0.85*min(N_components[old_labels[i]]),1.02*max(N_components[old_labels[i]]))
            # elif i == 0:
            #     # ax5[i//2][i%2].set_xlim(0.6*min(N_components[old_labels[i]]),1.02*max(N_components[old_labels[i]]))
            # else:
            #     # ax5[i//2][i%2].set_xlim(0.95*min(N_components[old_labels[i]]),1.02*max(N_components[old_labels[i]]))
            if i in [0,2]:
                ax5[i//2][i%2].set_ylabel(ylabel,fontsize = 22)
            else:
               ax5[i//2][i%2].yaxis.label.set_visible(False)
            ax5[i//2][i%2].set_xlabel(newlabels[i],fontsize = 22)            
            ax5[i//2][i%2].tick_params(axis = 'both', labelsize = 18)
    
    
        plt.subplots_adjust(wspace=0.13, hspace=0.15)
    
    # pid_df = pid_df.replace([np.inf, -np.inf], np.nan)
    # pid_df = pid_df.dropna(axis = 0)
    
    getFigure6666(pid_df, 'MAPE')
    
    
    
    
    
    
    ######Scatter plot of Model uncertainty vs Pearson r
    def getFigure4(data):
        fig2, ax2 = plt.subplots(3, 1, figsize=(11,24))
        data['total_error'] = data['H_in_situ']- data['I_crn_dca36']
        data['model_error'] = data['I_hv_temp_crn36'] - data['I_crn_dca36']
        data['random_error'] = data['total_error'] - data['model_error']
        
        data = data[['total_error', 'model_error', 'random_error', 
                     'NSE_36','color','landcover']]   
        min_total, max_total = min(data.total_error), max(data.total_error)
        min_model, max_model =  min(data.model_error), max(data.model_error)
        min_random, max_random =  min(data.random_error), max(data.random_error)
        
        R_mod = pearsonr(data['model_error'],data['NSE_36'])
        R_rnd = pearsonr(data['random_error'],data['NSE_36'])
        R_tot = pearsonr(data['total_error'],data['NSE_36'])
        Rs = [R_tot[0],R_mod[0] ,R_rnd[0]]
        print(Rs)
        print(data.head())
        sns.regplot(data['total_error'],data['NSE_36'],ax=ax2[0], 
                line_kws={"linestyle": '-.', 'color':'black'},
                )
        sns.regplot(data['model_error'],data['NSE_36'],ax=ax2[1], 
                line_kws={"linestyle": '-.', 'color':'black'},
                )
        sns.regplot(data['random_error'],data['NSE_36'],ax=ax2[2], 
                line_kws={"linestyle": '-.', 'color':'black'},
                )
        for i in data.values:
            print(i)
            ax2[0].scatter(i[0], i[3], c = i[-2], label = i[-1], s= 40)
            ax2[1].scatter(i[1], i[3], c = i[-2], label = i[-1], s= 40)
            ax2[2].scatter(i[2], i[3], c = i[-2], label = i[-1], s= 40)
        ax2[0].legend(*[*zip(*{l:h for h,l in zip(*ax2[0].get_legend_handles_labels())}.items())][::-1],
                      fontsize = 20, loc = 3)  
        ax2[1].legend(*[*zip(*{l:h for h,l in zip(*ax2[1].get_legend_handles_labels())}.items())][::-1],
                      fontsize = 20, loc = 3)  
        ax2[2].legend(*[*zip(*{l:h for h,l in zip(*ax2[2].get_legend_handles_labels())}.items())][::-1],
                      fontsize = 20, loc = 3)
        # [ax2[i].set_ylim(0, 1) for i in range(3)]
        ax2[0].set_xlim(0.995*min_total, max_total*1.005)
        ax2[1].set_xlim(0.96*min_model, max_model*1.007)   
        ax2[2].set_xlim(0.993*min_random, max_random*1.005)
        ax2[0].set_xlabel('SMAP informational total uncertainty, $I_{Tot}$',fontsize = 21)
        ax2[1].set_xlabel('SMAP informational model uncertainty, $I_{Mod}$',fontsize = 21)
        ax2[2].set_xlabel('SMAP informational random uncertainty, $I_{Rnd}$',fontsize = 21)
        [ax2[i].set_ylabel('NSE',fontsize = 21) 
          for i in range(3)]
        [ax2[i].tick_params(axis = 'both', labelsize = 18) for i in range(3)]
       
        [ax2[i].text(0.65,0.89, '$ r = {:.2f},\ p < {}$'.format(Rs[i],0.05), 
                  fontsize = 20,transform = ax2[i].transAxes) for i in np.arange(len(Rs))]
        plt.subplots_adjust(hspace = 0.16)
    getFigure4(pid_df)
    




















