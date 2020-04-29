#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 11:24:48 2020
@author: Ajit Johnson Nirmal
Adding metadata
"""

import pandas as pd


def add_meta (adata,meta,metafile_col,adata_col='ImageId'):
    """
    

    Parameters
    ----------
    adata : Ann Data Object
    meta : DataFrame
        CSV file containing the metadata information to be added to the adata object.
    metafile_col : string
        Column name in the metadata file that matches with one of the columns in adata.
    adata_col : string, optional
        Column name in adata object to match with the metadata file. The default is 'ImageId'.

    Returns
    -------
    adata : Modified AnnData Object
        DESCRIPTION.
        
    Example
    -------
    meta = pd.read_csv('/re_staining_tma/meta_data/meta_for_adata.csv')
    adata = add_meta (adata,meta,metafile_col='dearray_core_id',adata_col='ImageId')

    """
    
    # Find common elements between metafile_col and adata_col    
    a = list(meta[metafile_col])
    b = list(adata.obs[adata_col])        
    common_elements = list(set(a) & set(b))
    
    # Subset the meta file based on the common elements
    meta = meta[meta[metafile_col].isin(common_elements)]
    
    # adata OBS dataframe
    obs = pd.DataFrame(adata.obs[adata_col])
    obs['cellid'] = obs.index
    
    # Merge
    ob = obs.merge(meta, left_on=adata_col, right_on=metafile_col)
    ob.index = ob['cellid']
    ob = ob.reindex(obs.index)
    
    # Drop unwanted columns
    ob = ob.drop(['cellid', metafile_col, adata_col], axis=1)
    
    # Attach to adata
    for i in ob.columns:
        adata.obs[i] = ob[i]
        
    # Retun adata
    return adata

