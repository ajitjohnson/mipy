#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 12 09:18:09 2020
@author: Ajit Johnson Nirmal
local neighbourhood analysis
"""

import pandas as pd
import numpy as np
from sklearn.neighbors import BallTree
import multiprocessing as mp
import scipy
import ray 



def nn_interaction (adata,x_coordinate='X_position',y_coordinate='Y_position',phenotype='phenotype',
                    neighbour_radius=10, permutation=1000, p_val=0.05,image_id=None):
    
    if image_id is not None:
        adata = adata[adata.obs['ImageId'] == image_id] 
    
    # Create a dataFrame with the necessary inforamtion
    data = pd.DataFrame({'x': adata.obs[x_coordinate], 'y': adata.obs[y_coordinate], 'phenotype': adata.obs[phenotype]})
    
    # Identify neighbours
    print("Identifying neighbours within " + str(neighbour_radius) + " pixels of every cell")
    kdt = BallTree(data[['x','y']], metric='euclidean') 
    ind = kdt.query_radius(data[['x','y']], r=neighbour_radius)
    #remove self
    for i in range(0, len(ind)): ind[i] = np.delete(ind[i], np.argwhere(ind[i] == i))
    
    # Map phenotype
    phenomap = dict(zip(list(range(len(ind))), data['phenotype'])) # Used for mapping
    neighbours = pd.DataFrame(ind.tolist(), index = data.index) # neighbour DF
    
    # Loop through (all functionized methods were very slow)
    print("Mapping phenotype to neighbors")
    for i in neighbours.columns:
        neighbours[i] = neighbours[i].dropna().map(phenomap, na_action='ignore')
    
    # Drop NA
    neighbours = neighbours.dropna(how='all')
    
    # Collapse all the neighbours into a single column
    n = pd.DataFrame(neighbours.stack(), columns = ["neighbour_phenotype"])
    n.index = n.index.get_level_values(0) # Drop the multi index
    
    # Merge with real phenotype
    n = n.merge(data['phenotype'], how='inner', left_index=True, right_index=True)
    
    # Permutation
    print('Performing '+ str(permutation) + ' permutations')
    ray.init(ignore_reinit_error=True)
    @ray.remote
    def permutation_pval (data):
        data['neighbour_phenotype'] = np.random.permutation(data['neighbour_phenotype'])
        data_freq = data.groupby(['phenotype','neighbour_phenotype']).size().unstack()
        data_freq = data_freq.fillna(0).stack().values 
        return data_freq
    futures = [permutation_pval.remote(data=n) for i in range(permutation)]
    perm = pd.DataFrame((ray.get(futures))).T
    ray.shutdown()
    
    print('Consolidating tbe permutation results')
    # Calculate P value
    n_freq = n.groupby(['phenotype','neighbour_phenotype']).size().unstack().fillna(0).stack() 
    mean = perm.mean(axis=1)
    std = perm.std(axis=1)
    # Calculation
    z_scores = (n_freq.values - mean) / std
    z_scores[np.isnan(z_scores)] = 0
    p_values = scipy.stats.norm.sf(abs(z_scores))*2
    p_values = p_values[~np.isnan(p_values)]
    
    # Normalize the data by dividing each row by the max value
    #k = n.groupby(['phenotype','neighbour_phenotype']).size().unstack().fillna(0)
    #k_max = k.div(k.max(axis=1), axis=0).stack()
    
    # Normalize based on total cell count
    k = n.groupby(['phenotype','neighbour_phenotype']).size().unstack().fillna(0)
    total_cell_count = data['phenotype'].value_counts()
    total_cell_count = total_cell_count.reindex(k.columns).values    
    k_max = k.div(total_cell_count, axis = 0)
    k_max = k_max.div(k_max.max(axis=1), axis=0).stack()
    
    # DataFrame with the neighbour frequency and P values
    neighbours = pd.DataFrame({'count': k_max.values,'p_val': p_values}, index = k_max.index)
    neighbours.loc[neighbours[neighbours['p_val'] > p_val].index,'count'] = np.NaN
    del neighbours['p_val']
    neighbours = neighbours['count'].unstack()
        
    # Add to anndata
    adata.uns['local_neighbourhood'] = neighbours
    
    # return
    return adata

