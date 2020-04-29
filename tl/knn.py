#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 11:15:19 2020
@author: Ajit Johnson Nirmal
K Nearest Neighbour analysis
What cell type occurs most frequently close to the cell type of interest?
"""

import ray
import pandas as pd
import numpy as np
from sklearn.neighbors import BallTree
import scipy
#import setproctitle


def knn (adata,x_coordinate='X_position',y_coordinate='Y_position',phenotype='phenotype',k=2,
         permutation=1000,p_val=0.05,image_id=None):
    """
    

    Parameters
    ----------
    adata : AnnData Object
    x_coordinate : string, required
        Column name of the x-coordinates. The default is 'X_position'.
    y_coordinate : string, required
        Column name of the y-coordinates. The default is 'Y_position'.
    phenotype : string, required
        Column name of the column containing the phenotype information. The default is 'phenotype'.
    permutation : int, optional
        Number of permutations to calculate p-value. The default is 1000.
    p_val : float, optional
        The threshold below which will be considered as significant observations. The default is 0.05.
    image_id : string, optional
        For datasets containing multiple images, column name of the column containing the image id need to be passed here. The default is None.

    Returns
    -------
    adata : AnnData Object
        Returns the updated AnnData object with the results stored in adata.uns['knn'].
    
    Example
    -------
    adata = knn (adata,,phenotype='phenotype',permutation=1000,p_val=0.05)


    """
    
    print ('Identifying the Nearest Neighbours')
    # Create a data frame with the required data for this analysis
    data = pd.DataFrame({'x': adata.obs[x_coordinate], 
                         'y': adata.obs[y_coordinate], 
                         'phenotype': adata.obs[phenotype]})
    
    #Nearest neighbour identifier 
    tree = BallTree(data[['x','y']], leaf_size= 2)
    dist, ind = tree.query(data[['x','y']], k=k)
    
    # Map phenotype
    phenomap = dict(zip(list(range(len(ind))), data['phenotype'])) # Used for mapping
    neighbours = pd.DataFrame(ind.tolist(), index = data.index) # neighbour DF
    neighbours.drop(0, axis=1, inplace=True) # Remove self neighbour
    
    # Loop through (all functionized methods were very slow)
    print("Mapping phenotype to neighbors")
    for i in neighbours.columns:
        neighbours[i] = neighbours[i].dropna().map(phenomap, na_action='ignore')
    
    # Drop NA
    neighbours = neighbours.dropna(how='all')
    
    # Collapse all the neighbours into a single column
    knn = pd.DataFrame(neighbours.stack(), columns = ["neighbour_phenotype"])
    knn.index = knn.index.get_level_values(0) # Drop the multi index
    
    # Merge with real phenotype
    knn = knn.merge(data['phenotype'], how='inner', left_index=True, right_index=True)
    knn.columns = ['neighbour_phenotype','cell_phenotype']
    
    # Master neighbour dataFrame
    #knn = pd.DataFrame({'cell_phenotype': list(data['phenotype']),
    #                    'neighbour_phenotype': list(data.reset_index(drop=True).reindex(ind[:,1])['phenotype']),
    #                    'distance': list(dist[:,1])})
    
    #Freq Square Form one liner
    #knn_freq = knn.groupby(['cell_phenotype','neighbour_phenotype']).size().unstack()
    #knn_distance = knn.groupby(['cell_phenotype','neighbour_phenotype'])['distance'].median().unstack()
    
    # Permutation
    print('Performing '+ str(permutation) + ' permutations')
    ray.init(ignore_reinit_error=True)
    @ray.remote
    def permutation_pval (data):
        data['neighbour_phenotype'] = np.random.permutation(data['neighbour_phenotype'])
        data_freq = data.groupby(['cell_phenotype','neighbour_phenotype']).size().unstack()
        data_freq = data_freq.fillna(0).stack().values 
        return data_freq
    futures = [permutation_pval.remote(data=knn) for i in range(permutation)]
    perm = pd.DataFrame((ray.get(futures))).T
    ray.shutdown()
    
    print('Consolidating tbe permutation results')
    # Calculate P value
    knn_freq = knn.groupby(['cell_phenotype','neighbour_phenotype']).size().unstack().fillna(0).stack() 
    mean = perm.mean(axis=1)
    std = perm.std(axis=1)
    # Calculation
    z_scores = (knn_freq.values - mean) / std
    p_values = scipy.stats.norm.sf(abs(z_scores))*2
    p_values = p_values[~np.isnan(p_values)]
        
    # Normalize the data by dividing each row by the max value
    #k = knn.groupby(['cell_phenotype','neighbour_phenotype']).size().unstack().fillna(0)
    #k_max = k.div(k.max(axis=1), axis=0).stack()
    
    # Normalize based on total cell count
    k = knn.groupby(['cell_phenotype','neighbour_phenotype']).size().unstack().fillna(0)
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
    adata.uns['knn'] = neighbours
    
    # return
    return adata