# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 15:15:03 2020
@author: Ajit Johnson Nirmal
Plot to vizualise the effect of background normalization
"""

# Import library
import seaborn as sns; sns.set(color_codes=True)
import matplotlib.pyplot as plt
import numpy as np


# Function
def raw_vs_normalized_plot (adata, out_dir, log=True):
    '''
    Parameters:
        adata- AnnData object created using the function df_to_annobject.
        out_dir: Directory to which the images should be saved.
        log: Bool. If True, then the data will be converted to log scale.
    Returns:
        PDF plots that will be saved in the specified directory
    Example:
        raw_vs_normalized_plot (adata, out_dir = "C:/Users/ajit/plots", log=True)
    '''
    # Figure
    def plot_normalization_figure (adata, marker, out_dir):
        # Data for plotting
        m_idx = adata.var.index.tolist().index(marker) # Get the index of marker of interest
        # Raw data
        data = np.array(adata.raw.X[m_idx]).flatten()
        # Normalized data 
        n_data = np.array(adata.X[:,m_idx]).flatten()
        # Plot
        sns.set_style("white")
        if log == True:
            sns.distplot(np.log1p(data), hist=False, rug=False, label="Before Normalization");     
            sns.distplot(np.log1p(n_data), hist=False, rug=False, label="After Normalization");
        else:
            sns.distplot(data, hist=False, rug=False, label="Before Normalization");     
            sns.distplot(n_data, hist=False, rug=False, label="After Normalization");
        plt.savefig(out_dir + "/" + marker + ".pdf")
        plt.clf()
    # Run the function
    r_figure = lambda x: plot_normalization_figure(marker = x, adata = adata, out_dir=out_dir)
    figures = list(map(r_figure, adata.var.index)) # Apply function
        
