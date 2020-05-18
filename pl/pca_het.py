#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 09:28:00 2020
@author: Ajit Johnson Nirmal
PCA heterogenity plot
"""

from sklearn.decomposition import PCA
from sklearn.preprocessing import MinMaxScaler
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import seaborn as sns; sns.set(color_codes=True)

def pca_het (adata, x_coordinate='X_position',y_coordinate='Y_position',
             unique_id='ImageId',subset_image=None, raw_data=True,
             phenotype='phenotype',phenotype_of_interest=None,
             genes_of_interest=None,s = 2, alpha=0.8,fontsize=5,
             widths = [2],heights = [6,1], save_fig=False, save_dir=None,
             save_name = 'PCA_Heterogenity_Plot', save_format='png', figsize=(10, 10)):
    
    
    # Copy adata
    bdata = adata
    
    if bdata.raw is None:
        bdata.raw = bdata

    if subset_image is not None:
        bdata = bdata[bdata.obs[unique_id] == subset_image]
    
    if phenotype_of_interest is not None:
        bdata = bdata[bdata.obs[phenotype].isin(phenotype_of_interest)]
    
    if genes_of_interest is not None:
        bdata = bdata[:, genes_of_interest]
    else:
        genes_of_interest = list(bdata.var.index)
    
           
    # Create a dataframe with the necessary information
    if raw_data is True:
        data = pd.DataFrame(np.log1p(bdata.raw[:, genes_of_interest].X), index=bdata.obs.index, columns=bdata.var.index)
    else:
        data = pd.DataFrame(bdata[:, genes_of_interest].X, index=bdata.obs.index, columns=bdata.var.index)
    
    coordinates = pd.DataFrame({'x':bdata.obs[x_coordinate].values, 'y':bdata.obs[y_coordinate].values},index=bdata.obs.index)
    
    pca = PCA(n_components=3)
    pca.fit(data)
    
    # For heatmap
    ev_m = pd.DataFrame(pca.components_, columns= data.columns, index=['R','G','B'])

    # For scatter plot
    ev_cell = pd.DataFrame(pca.transform(data))
    ev_cell.index = data.index
    ev_cell[0]= ev_cell[0].clip(np.percentile(ev_cell[0], 2.5),np.percentile(ev_cell[0], 97.5))
    ev_cell[1]= ev_cell[1].clip(np.percentile(ev_cell[1], 2.5),np.percentile(ev_cell[1], 97.5))
    ev_cell[2]= ev_cell[2].clip(np.percentile(ev_cell[2], 2.5),np.percentile(ev_cell[2], 97.5))
    
    # Coloring
    scaler = MinMaxScaler(feature_range=(0, 255))
    scaler.fit(ev_cell)
    ev_cell_rgb = scaler.transform(ev_cell)
    ev_cell_rgb = ev_cell_rgb.astype(int)
    
    # Plot 
    fig = plt.figure(constrained_layout=True,figsize=figsize)
    spec = gridspec.GridSpec(ncols=1, nrows=2, figure=fig, width_ratios=widths,height_ratios=heights)
    ax1 = fig.add_subplot(spec[0, 0])
    ax2 = fig.add_subplot(spec[1, 0])
    ax1.set_facecolor('black')
    ax1.scatter(adata.obs[x_coordinate].values, 
                adata.obs[y_coordinate].values, 
                c='#282828', s=s, alpha=1, marker=",")
    ax1.scatter(bdata.obs[x_coordinate].values, 
                bdata.obs[y_coordinate].values, 
                c=ev_cell_rgb/255, s=s, alpha=alpha,marker=",")
    ax1.axes.get_xaxis().set_ticks([])
    ax1.axes.get_yaxis().set_ticks([])
    ax1.invert_yaxis()
    heatmap = sns.heatmap(ev_m, cmap='vlag',ax=ax2, xticklabels = 1, yticklabels = 1) #standard_scale=1 , z_score=0
    heatmap.set_xticklabels(heatmap.get_xticklabels(), rotation=90,fontweight='bold',fontsize=fontsize)
    heatmap.set_yticklabels(heatmap.get_yticklabels(), fontweight='bold',fontsize=fontsize)
    plt.yticks(rotation=0) 
    
    if save_fig is True:
        if save_dir is not None:
            plt.savefig(save_dir + "/" + save_name + "." + save_format, bbox_inches="tight", dpi=300)
        else:
            plt.savefig(save_name + "." + save_format, bbox_inches="tight", dpi=300)
    else:
        plt.show()
            
        


