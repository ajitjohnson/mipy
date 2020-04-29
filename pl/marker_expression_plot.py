# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 11:23:03 2020
@author: Ajit Johnson Nirmal
Plots of expression vs prediction
"""

# Import library
import matplotlib.pyplot as plt
import seaborn as sns; sns.set(color_codes=True)
import numpy as np

def marker_expression_plot (adata, x='X_position', y='Y_position', gate = 0.5, color=None, save_figure= False, plot_all=False, out_dir=None, alpha= 0.3, figure_type="png", s=1, cmap='gist_heat', figsize=(12, 5)):
    '''
    Parameters:
        adata- AnnData object containing both expression and marker expression predictions.
        color- Marker to color the cells by.
        x- column name of x-coordinates.
        y- column name of y-coordinates.
        s- size of the plot points.
        gate: int. Default is set to 0.5. By default rescale function, scales the data such that values above 0.5 are considered positive cells.
        out_dir: Directory to which the figures will be saved. If no directory is passed, figures will be saved in the working directory.
        figure_type: Png or Pdf.
        plot_all: Plot all the markers and save it the given directory.
    Returns:
        scatter plot
    Example:
        plots = marker_expression_plot (adata, x='X_position', y='Y_position', 
                                    color=None, save_figure= True, plot_all=True, 
                                    out_dir=out_dir, alpha= 0.3, figure_type="png", 
                                    s=1, cmap='gist_heat', figsize=(12, 5))
    '''  
    # Figure settings
    sns.set_style("white")
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    
    x= adata.obs[x]
    y= adata.obs[y]
    
    if plot_all == True:
        for i in adata.var.index:
            # Print on screen
            print('Plotting ' + str(i))
            # Identify the marker to plot
            c = adata[:,i].X
            c_col = ['red' if x >= gate else 'black' for x in c]
            c_exp = np.concatenate(adata.raw[:,i].X, axis=0 )
            
            # Plot all figures in a loop
            fig, axes = plt.subplots(nrows=1, ncols=2,figsize=figsize)
            axes[0].scatter(x, y, c= np.log1p(c_exp), edgecolors='none', s=s, cmap=cmap, alpha=alpha)
            plt.title(i)
            axes[0].invert_yaxis()
            axes[1].scatter(x, y, c=c_col, edgecolors='none', s=s, alpha=alpha)
            plt.title(i)
            axes[1].invert_yaxis()
            fig.tight_layout()
            if save_figure == True:
                if out_dir != None:
                    plt.savefig(out_dir + "/" + i + "." + figure_type, bbox_inches="tight")
                else:
                    plt.savefig(i + "." + figure_type, bbox_inches="tight")
            else:
                plt.show()
            # Clear plot
            plt.clf()
    else:
        # Plot a single figure
        c = adata[:,color].X
        c_col = ['red' if x >= gate else 'black' for x in c]
        c_exp = np.concatenate(adata.raw[:,color].X, axis=0 )
        
        # Plot
        fig, axes = plt.subplots(nrows=1, ncols=2,figsize=figsize)
        axes[0].scatter(x, y, c= np.log1p(c_exp), edgecolors='none', s=s, cmap=cmap, alpha=alpha)
        plt.title(color)
        axes[0].invert_yaxis()
        axes[1].scatter(x, y, c=c_col, edgecolors='none', s=s, alpha=alpha)
        plt.title(color)
        axes[1].invert_yaxis()
        fig.tight_layout()
        if save_figure == True:
            if out_dir != None:
                plt.savefig(out_dir + "/" + color + "." + figure_type, bbox_inches="tight")
            else:
                plt.savefig(color + "." + figure_type, bbox_inches="tight")
        else:
            plt.show()
        # Clear plot
        plt.clf()
            

    


            
        
        
        