#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 11 13:52:37 2020
@author: Ajit Johnson Nirmal
Merge list of Ann Data objects
"""

import scanpy as sc

def adata_concatenate (adata_loc_list):
    """
    

    Parameters
    ----------
    adata_loc_list : list
        A list containing the location to the Ann Data objects to be concatenated.

    Returns
    -------
    adata : Ann Data Object
        Concatenated Ann Data.
    
    Example
    -------
    adata_loc_list = ['/Users/aj/data/adata1.h5ad',
                      '/Users/aj/data/adata2.h5ad']
    adata = adata_concatenate (adata_loc_list)


    """
    
    def loc_to_adata (data):
        # Print
        print('Reading: ' + data)
        # read adata
        adata = sc.read(data)
        return adata
    
    # Apply function
    r_loc_to_adata = lambda x: loc_to_adata(data=x) # Create lamda function 
    all_adata = list(map(r_loc_to_adata, adata_loc_list)) # Apply function 
    
    # Concatenate the AnnData objects
    adata = all_adata[0].concatenate(all_adata[1:], join='outer',index_unique=None)
    # Add the uns elements
    if all_adata[0].uns is not None:
        adata.uns = all_adata[0].uns
    
    
    # return results
    return adata