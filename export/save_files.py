#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 19:28:38 2020
@author: Ajit Johnson Nirmal
Functions to export results
"""

import pandas as pd

def save_phenotype (adata,phenotype_column,file_name=None,out_dir=None):
    """
    

    Parameters
    ----------
    adata : AnnData Object
    phenotype_column : string
        Name of the column that contains the phenotype information.
    file_name : String, optional
        Name the file. The default is None.
    out_dir : string, optional
        Location to save the CSV file. The default is Working Directory.

    Returns
    -------
    Exports a CSV file.
    
    Example
    -------
    save_phenotype (adata,phenotype_column='kmeans_renamed',
    file_name='PTCL phenotype',out_dir=''/Users/aj/PTCL1')

    """
    
    # Create a dataframe 
    save_data = adata.obs[phenotype_column]
    
    # Add File Name
    if file_name==None:
        file_name = phenotype_column
    
    #Save as CSV
    if out_dir == None:
        save_data.to_csv(file_name + '.csv')
    else:
        save_data.to_csv(out_dir + '/' + file_name + '.csv')
        

def save_phenotype_proportion (adata,phenotype_column,file_name=None,out_dir=None):
    """
    
    Parameters
    ----------
    adata : AnnData Object
    phenotype_column : string
        Name of the column that contains the phenotype information.
    file_name : String, optional
        Name the file. The default is None.
    out_dir : string, optional
        Location to save the CSV file. The default is Working Directory.

    Returns
    -------
    Exports a CSV file.
    
    Example
    -------
    save_phenotype_proportion (adata,phenotype_column='kmeans_renamed',
    file_name='PTCL Proportion',out_dir=''/Users/aj/PTCL1')

    """
    
    # Create DataFrame
    prop_data = pd.DataFrame(adata.obs[phenotype_column].value_counts())
    prop_data['proportion_cells'] = round((prop_data[phenotype_column] / sum(prop_data[phenotype_column])) * 100, 2)
    prop_data.columns = ['total_cells','proportion_cells']
    
    # Add File Name
    if file_name==None:
        file_name = phenotype_column + '_proportion'
    
    #Save as CSV
    if out_dir == None:
        prop_data.to_csv(file_name + '.csv')
    else:
        prop_data.to_csv(out_dir + '/' + file_name + '.csv')
