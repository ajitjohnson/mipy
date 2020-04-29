# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 10:48:58 2020
@author: Ajit Johnson Nirmal
Identify mutually exclusive markers for background removal
"""

# Import library
import itertools
from  more_itertools import unique_everseen
from sklearn.utils.extmath import randomized_svd
import numpy as np
import pandas as pd
import scanpy as sc

def mutually_exclusive_markers (adata, top_features=3, log=True, sample_data=True, fraction=None, random_sample=None):
    '''
    Parameters:
        adata: AnnData object created using the function df_to_annobject
        top_features: int. Number of mutually exclusive marker to be identified. Defaults to 3.
        log: Boolian. Set to True to log the data before performing SVD. Defaults to True.
        sample_data: Boolian. If true, either use 'fraction' or 'random_sample' to randomly choose a subset of the dataset.
        fraction: use valuses between [0-1]. Subsample to this fraction of the number of observations.
        random_sample: int. Subsample to this number of observations.
    Returns:
        AnnData with the mutually exclusive markers: adata.uns['mutually_exclusive_markers']
        The SVD scores of the top mutually exclusive markers: adata.uns['mutually_exclusive_markers_with_scores']
        SVD scores of all the marker pairs: adata.uns['svd_ratio_matrix']: Use for generating heatmap.
    Example:
        mutually_exclusive_markers (adata, top_features=5, log=True)
    '''
    # Find all pairs of markers
    markers= adata.var.index
    def combination_tuple (d, markers):
        res = [(d, val) for val in markers]
        return res
    # Apply function
    r_tuple = lambda x: combination_tuple(d=x, markers=markers)
    all_combinations = list(map(r_tuple, markers))
    all_combinations = list(itertools.chain.from_iterable(all_combinations))
   # all_combinations = list(combinations_with_replacement(adata.var.index,2))
   
    # Randomly sample data
    if sample_data == True:
        bdata = sc.pp.subsample(adata,fraction=fraction, n_obs=random_sample,copy=True)
    else:
        bdata = adata
    
    # Function to perform SVD on a single pair of markers
    def svd (d, comb, log):
        # Subset the data based on the particular marker combination    
        data_sub = [d[:,comb[0]].X, d[:,comb[1]].X]
        if log == True:
            data_sub = np.log10(data_sub)
        u, s, vt = randomized_svd(np.array(data_sub), n_components = 2)
        r = s[1]/s[0]
        return r # return the ratio
    
    # Run the SVD function on all marker pairs (vectorized)
    r_lamda = lambda x: svd(comb=x, d=bdata, log=log) # Create lamda function 
    all_r = list(map(r_lamda, all_combinations)) # Apply function
    
    # Create a dataframe with the results
    exclusive_markers = pd.DataFrame(all_combinations)
    exclusive_markers.columns = ['MarkerA', 'MarkerB']
    exclusive_markers['r'] = all_r
    
    # Spare matrix of the SVD ratio
    #svd_ratio = pd.DataFrame(ssd.squareform(all_r))
    svd_ratio = pd.DataFrame(np.array_split(np.array(all_r),len(markers)))
    svd_ratio.index = list(unique_everseen(exclusive_markers.iloc[:,0]))
    svd_ratio.columns = list(unique_everseen(exclusive_markers.iloc[:,0]))
    
    # Sort the markers based on the Marker A column and the r values and get the top hits
    df1 = exclusive_markers.groupby(["MarkerA"])    
    df2 = df1.apply (lambda x: x.sort_values(['r'], ascending=False))
    df3=df2.reset_index(drop=True)
    ex_markers = df3.groupby('MarkerA').head(top_features)
    
    # Add the identified markers as a seperate column
    ex_markers_table = ex_markers.groupby('MarkerA')['MarkerB'].apply(' '.join).reset_index()
    table_elements = ex_markers_table['MarkerB'].apply(lambda x: pd.Series(x.split(' ')))
    ex_markers_table = pd.DataFrame(ex_markers_table['MarkerA']).join(table_elements)  
    
    # Add the elements to AnnData
    adata.uns['mutually_exclusive_markers'] = ex_markers_table
    adata.uns['mutually_exclusive_markers_with_scores'] = ex_markers
    adata.uns['svd_ratio_matrix'] = svd_ratio
    
    # Return the final dataframe
    return adata