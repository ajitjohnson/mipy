#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 19:51:38 2020
@author: Ajit Johnson Nirmal
Clustering and Visualisation of nn_communities
"""

import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np


def community_clustering (adata, community_matrix='nn_communities',subset_image=None,unique_id='ImageId',
                          n_neighbors=30, method='leiden',resolution=1,random_state=0,label=None,smallest_cluster=None):
    """
    

    Parameters
    ----------
    adata : Ann Data Object
    community_matrix : DataFrame, optional
        Resultant DataFrame of running the function nn_communities. The default is 'nn_communities'.
    subset_image : string, optional
        ImageId of a single image to be subsetted for analyis. The default is None.
    unique_id : string, optional
        For datasets containing multiple images, column name of the column containing the image id 
        need to be passed here. The default is 'ImageId'.
    n_neighbors : int, optional
        The size of local neighborhood (in terms of number of neighboring data points) 
        used for manifold approximation. Larger values result in more global views of the 
        manifold, while smaller values result in more local data being preserved. In general 
        values should be in the range 2 to 100. The default is 30.
    method : string, optional
        Available option a) leiden b) louvain c) phenograph. The default is 'leiden'.
    resolution : int, optional
        A parameter value controlling the coarseness of the clustering. Higher values lead to more clusters. The default is 1.
    random_state : int, optional
        Change the initialization of the optimization. The default is 0.
    label : string, optional
        Name the resulting column that is returned. The default is the name of the method used.
    smallest_cluster : int, optional
        Merge small clusters into a single cluster. Any cluster containing fewer cells than 
        the provided number will be grouped into a single category called 'trimmed_clusters'. The default is None.

    Returns
    -------
    adata : Modified Ann Data Object
        Check adata.obs.
        
    Example
    -------
    adata = nn_clustering (adata, unique_id='ImageId',n_neighbors=30, 
                           method='leiden',resolution=0.5)

    """
    
    # Remove previous iteration of clustering
    if label is not None:
        if label in adata.obs.columns:
            del adata.obs[label]
            del adata.uns[str(label)+'_colors']
    elif method in adata.obs.columns:
        del adata.obs[method]
        del adata.uns[str(method)+'_colors']
        
    
    # Subset a particular image if needed
    if subset_image is not None:
        bdata = adata[adata.obs[unique_id] == subset_image]
    else:
        bdata = adata
    
    d = bdata.uns[community_matrix].fillna(0)   
    d = d[(d.T != 0).any()] # Drop all rows with zero   
    samples = d.index  # Get the cells that survived
    cdata = bdata[samples] # Do this so that the .obs of this can be used below
    
    # Create a new ann data object
    bdata = ad.AnnData(d)
    bdata.obs = pd.DataFrame(cdata.obs)
    
    # Cluster data
    sc.tl.pca(bdata)
    sc.pp.neighbors(bdata, n_neighbors=n_neighbors)
    
    # Method of clustering
    if method == 'leiden':
        sc.tl.leiden(bdata,resolution=resolution,random_state=random_state)
    if method == 'louvain':
        sc.tl.louvain(bdata,resolution=resolution,random_state=random_state)
        
    # Collapse small clusters into a single cluster
    if smallest_cluster is not None:
        passed_clusters = list(bdata.obs[method].value_counts()[bdata.obs[method].value_counts() < smallest_cluster].index)
        # create dict
        mapping = dict(zip(passed_clusters, np.repeat('trimmed_clusters',len(passed_clusters))))
        bdata.obs[method] = bdata.obs[method].replace(mapping)
        bdata.obs[method] = bdata.obs[method].astype('category')
        
        
    # Merge the clustering results with adata
    original = adata.obs
    results = bdata.obs[method]
    r =  original.merge(results, how='outer', left_index=True, right_index=True)
    r = r.reindex(adata.obs.index) # Reindex
    
    # Modify adata
    if label is not None:
        adata.obs[label] = r[method]
    else:
        adata.obs[method] = r[method]
        
    # Return
    return adata
        
    
