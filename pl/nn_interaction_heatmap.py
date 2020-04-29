#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 16:34:23 2020
@author: Ajit Johnson Nirmal
Neighbourhood Analysis HeatMap
"""

from scipy.cluster.hierarchy import dendrogram, linkage
import seaborn as sns



def nn_interaction_heatmap (adata, neighbourhood_result):
    """
    

    Parameters
    ----------
    adata : Ann Data Object
    neighbourhood_result : string
        Name of neighbourhood function that was performed.

    Returns
    -------
    HeatMap.
    
    Example
    -------
    neighbourhood_heatmap (adata, neighbourhood_result='local_neighbourhood')

    """
    
    
    neighbours_max = adata.uns[neighbourhood_result]
    neighbours_max = neighbours_max.fillna(0)
    
    # Cluster the data for better Viz
    Z = linkage(neighbours_max, 'ward')   
    clust = dendrogram(Z, no_plot=True).get("ivl")
    
    # Re format the neighbours_max matrix with the clustrring
    index_name = []
    for i in clust:
        index_name.append(neighbours_max.iloc[int(i)].name)
    
    neighbours_max = neighbours_max.reindex(index_name)
    neighbours_max = neighbours_max[index_name]
    
    # Plot
    mask = neighbours_max.isnull()
    sns.set_style("whitegrid", {'axes.grid' : True})
    sns.clustermap(neighbours_max, cmap = 'Blues', mask=mask, row_cluster=False, col_cluster=False)
    






