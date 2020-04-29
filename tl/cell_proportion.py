#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 21:37:31 2020
@author: Ajit Johnson Nirmal
Regional Cell Proportion Analysis
"""

import pandas as pd

def cell_proportion (adata, regions, phenotype):
    
    p = pd.DataFrame()
    for r in regions:
        prop = pd.DataFrame(adata[adata.obs[r] == r].obs[phenotype].value_counts())
        prop.columns = [r]    
        p = p.merge(prop, how='outer', left_index=True, right_index=True)
    
    # impute 0 to nan
    proportion = p.fillna(0)
    
    # Convert to proportion
    #proportion = proportion / proportion.sum(axis = 0)
    
    
    # Attach to adata
    adata.uns['cell_proportion'] = proportion
    
    # Return
    return adata