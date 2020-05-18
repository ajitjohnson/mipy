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

def pca_het (adata, x_coordinate='X_position',y_coordinate='Y_position',
             unique_id='ImageId',subset_image=None,
             phenotype='phenotype',phenotype_of_interest=None,
             genes_of_interest=None):
    
    # Copy adata
    bdata = adata
    
    if subset_image is not None:
        bdata = bdata[bdata.obs[unique_id] == subset_image]
    
    if phenotype_of_interest is not None:
        bdata = bdata[bdata.obs[phenotype].isin(phenotype_of_interest)]
    
    if genes_of_interest is not None:
        bdata = bdata[:, genes_of_interest]
        
        
    # Create a dataframe with the necessary information
    data = pd.DataFrame(np.log1p(bdata.raw.X), index=bdata.obs.index, columns=bdata.var.index)
    coordinates = pd.DataFrame({'x':bdata.obs[x_coordinate].values, 'y':bdata.obs[y_coordinate].values},index=bdata.obs.index)
    
    pca = PCA(n_components=3)
    pca.fit(data)
    ev_m = pd.DataFrame(pca.fit_transform(data.T))
    ev_m.index = data.columns
    
    ev_cell = pd.DataFrame(pca.fit_transform(data))
    ev_cell.index = data.index
    
    # Coloring
    scaler = MinMaxScaler(feature_range=(0, 255))
    scaler.fit(ev_cell)
    ev_cell_rgb = scaler.transform(ev_cell)
    ev_cell_rgb = ev_cell_rgb.astype(int)
    
    # Plot
    fig, ax = plt.subplots()
    plt.scatter(bdata.obs[x_coordinate].values, 
                bdata.obs[y_coordinate].values, 
                c=ev_cell_rgb/255.0, s = 2, alpha=0.8)
    ax.invert_yaxis() 









