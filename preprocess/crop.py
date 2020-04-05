#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 10:27:16 2020
@author: Ajit Johnson Nirmal
Cropping regions of interest
"""


def crop (adata):










a = np.where(adata.obs['kmeans_renamed'] == 'tumor- PD1+', 1, 0)
c_col = ['red' if x==1 else 'black' for x in a]








      
        marker = 'CD45'
        x = adata.obs['X_position']
        y = adata.obs['Y_position']
        
        m_idx = adata.var.index.tolist().index(marker) # Get the index of marker of interest
        tmp_dataframe = pd.DataFrame(adata.X)
        hue = np.array(tmp_dataframe[m_idx])
        hue = ['red' if x > 0.5 else 'black' for x in hue]
        
        # Plotting    
        fig, ax = plt.subplots()
        pts = ax.scatter(x, y, s=1,c=hue,cmap='gist_heat')
        ax.invert_yaxis()
        
        # Function call to do the Lasso selection
        selector = SelectFromCollection(ax, pts)
        # Return indeces of the selected points
        tumor_idx= selector.ind
        len(tumor_idx)
        
        # Update adata
        adata.obs.loc[adata.obs.index[tumor_idx], 'roi-1'] = "roi-1"
        adata.obs['roi-1'].value_counts() # Checking