#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 13:58:01 2020
@author: Ajit Johnson Nirmal
Infiltration/Proximity score for any given two groups
"""

import pandas as pd
import itertools

def proximal_density_score (adata, goi, community_matrix='nn_communities_raw',
                            label= 'pds'):
    """
    

    Parameters
    ----------
    adata : Ann Data Object
        DESCRIPTION.
    goi : list
        List of phrnotypes in community_matrix to be considered as a group.
    community_matrix : DataFrame, optional
        Resultant DataFrame of running the function nn_communities. The default is 'nn_communities'.
    label : string, optional
        Name the resulting column that is returned. The default is 'pds'.

    Returns
    -------
    adata : Modified Ann Data Object
        Check adata.obs['pds']
    
    Example
    -------
    goi = ['Tumor CD30+']
    adata = proximal_density_score (adata, goi, label= 'tumor_macs')

    """
    
    
    # Create the dataframe
    d = adata.uns[community_matrix].fillna(0)
    
    # Subset the groups of interest
    def pds (d, comb, goi):
        gA= pd.DataFrame(pd.DataFrame(d[goi]).sum(axis=1))
        gB= pd.DataFrame(pd.DataFrame(d[comb]).sum(axis=1))
        d = gA.merge(gB,how='outer', left_index=True, right_index=True)
        d = d.div(d.sum(axis=1), axis=0).fillna(0) 
        d.iloc[:, 0] = d.iloc[:, 0] * -1
        d[str(comb)]  = d.sum(axis=1)
        d = pd.DataFrame(d.iloc[:, 2])
        return d
    
    all_comb = d. columns
    # Run function
    r_pds = lambda x: pds(d=d, goi=goi, comb=x) # Create lamda function 
    all_data = list(map(r_pds, all_comb)) # Apply function 
    
    result = pd.concat(all_data, join='outer', axis=1) 
    
    # Reindex
    result = result.reindex(adata.obs.index)

    # Return
    adata.uns[label] = result
    
    return adata
