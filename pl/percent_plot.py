#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 14:00:04 2020
@author: Ajit Johnson Nirmal
Stacked Bar Plot for visualizing proportion of cells
"""

import pandas as pd
import matplotlib.pylab as plt
import seaborn as sns; sns.set(color_codes=True)
sns.set(style="white")

def percent_plot (adata,x_axis,y_axis,method='percent',subset_data=None, figsize=(10, 10)):
    """
    

    Parameters
    ----------
    adata : AnnData object
    x_axis : string
        Column in adata that need to be plotted in the x-axis.
    y_axis : string
        Column in adata that need to be plotted in the y-axis.
    method : string, optional
        Available options: 'percent' and 'absolute'. Use Percent to plot as percent proportion.
        Use 'absolute' to plot the plot the absolute number.
        The default is 'percent'.
    subset_data: dict, optional
        Subset the data before plotting the data, where the dict key is the interested column name in adata and
        the values are the sample groups to be subsetted.
        The default is None
    figsize : tuple, optional
        Pass in two values for figure height and width. The default is (10, 10).

    Returns
    -------
    None.
    
    Example
    -------
    percent_plot (adata,x_axis='ImageId',y_axis='phenotype',method='percent',figsize=(10, 10))

    """
    
    data = pd.DataFrame(adata.obs)
    
    if subset_data is not None:     
        data = data[data[list(subset_data.keys())[0]].isin(list(subset_data.values())[0])]
        
    # Make a dataframe with the plot information
    r = data[[x_axis,y_axis]]    
    r[y_axis] = r[y_axis].astype('str')
    r[x_axis] = r[x_axis].astype('str')
    
    
    # Method: Absolute or Percentile
    if method == 'percent':
        total = r.groupby([x_axis,y_axis]).size().unstack().fillna(0).sum(axis=1)
        rg = pd.DataFrame(r.groupby([x_axis,y_axis]).size().unstack().fillna(0).div(total, axis=0).stack())
    else:
        rg = pd.DataFrame(r.groupby([x_axis,y_axis]).size().unstack().fillna(0).stack())

    # Change column name
    rg.columns = ['count']
    
    # Add the index as columns in the data frame    
    rg[x_axis] = rg.index.get_level_values(x_axis)
    rg[y_axis] = rg.index.get_level_values(y_axis)
    pivot_df = rg.pivot(index=x_axis, columns=y_axis, values='count')
    
    # Create color panel
    if len(rg[y_axis].unique()) <= 9:
        cmap = "Set1"        
    elif len(rg[y_axis].unique()) > 9 and len(rg[y_axis].unique()) <=20:
        cmap = plt.cm.tab20      #tab20  
    else:
        cmap = plt.cm.gist_ncar
        
    
    # Plotting
    p = pivot_df.plot.bar(stacked=True, cmap=cmap,figsize=figsize,width=.9)
    p.legend(bbox_to_anchor=(1.04,1), loc="upper left")
    
    
    
    
