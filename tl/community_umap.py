#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 23:24:38 2020
@author: Ajit Johnson Nirmal
UMAP of the cell communities
"""

import anndata as ad
import scanpy as sc
import pandas as pd


def community_umap (adata,community_matrix='nn_communities', n_neighbors=30, random_state=0):
    """
    

    Parameters
    ----------
    adata : Ann Data Object
    community_matrix : TYPE, optional
        Resultant DataFrame of running the function nn_communities. The default is 'nn_communities'.
    n_neighbors : int, optional
        The size of local neighborhood (in terms of number of neighboring data points) 
        used for manifold approximation. Larger values result in more global views of the 
        manifold, while smaller values result in more local data being preserved. In general 
        values should be in the range 2 to 100. The default is 30.
    random_state : int, optional
        Change the initialization of the optimization. The default is 0.

    Returns
    -------
    adata : Modified Ann Data Object
        Check adata.obsm['X_umap'].
    
    Example
    -------
    adata = community_umap (adata,community_matrix='nn_communities', n_neighbors=50)
    

    """
    
    # Create dataframe
    d = adata.uns[community_matrix].fillna(0) 
    
    # Create a new ann data object
    bdata = ad.AnnData(d)
        
    # Build the UMAP
    sc.tl.pca(bdata)
    sc.pp.neighbors(bdata, n_neighbors=n_neighbors)
    sc.tl.umap(bdata,random_state=random_state)
    
    # Add the UMAP embedding to adata
    adata.obsm['X_umap'] = bdata.obsm['X_umap']
    
    # Return
    return adata