#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 21:57:54 2020
@author: Ajit Johnson Nirmal
Using Napari to Visualize images overlayed with phenotypes
"""

#%gui qt
#from skimage import io
import napari
import pandas as pd
import random
#from skimage.external.tifffile import imread
#import numpy as np
import tifffile as tiff



def view_image (image_path, adata, phenotype_column='phenotype', 
                overlay_phenotype=None,markers=None,
                channel_names='default',
                x='X_position',y='Y_position',point_size=10,
                point_color=None):
    """
    Parameters
    ----------
    image_path : string
        Location to the image file.
    adata : AnnData Object
    phenotype_column : string, optional
        Name of the column that contains the phenotype information. The default is 'phenotype'.
    overlay_phenotype : list, optional
        If only specfic phenotypes need to be overlayed on the image, pass their names as a list. The default is None.
    markers : TYPE, optional
        DESCRIPTION. The default is None.
    channel_names : list, optional
        List of channels that need to be included. The default is adata.uns['all_markers'].
    x : TYPE, optional
        DESCRIPTION. The default is 'X_position'.
    y : TYPE, optional
        DESCRIPTION. The default is 'Y_position'.
    point_size : TYPE, optional
        DESCRIPTION. The default is 5.

    Returns
    -------
    None.

    """
    
    # Recover the channel names from adata
    if channel_names is 'default':
        channel_names = adata.uns['all_markers']
    else:
        channel_names = channel_names
    
    
    # Index of the marker of interest and corresponding names
    if markers is None:
        idx = list(range(len(channel_names)))
        channel_names = channel_names
    else:
        idx = []  
        for i in markers:
            idx.append(list(channel_names).index(i))
        channel_names = markers
    
        
    # Load the image
    #image = imread(image_path, key = idx)
    #image = io.imread(image_path, key = idx).T
    image = tiff.imread(image_path, key = idx)

    
    #def load_image (image_path, key):
    #    image = imread(image_path, key=key)
    #    return image
    #r_load_image = lambda x: load_image(image_path=image_path, key=x) # Create lamda function    
    #results = list(map(r_load_image, idx)) # Apply function 
    #image  = np.array(results)

        
    # Load the viewer
    viewer = napari.view_image(
    image,
    is_pyramid=False,
    channel_axis=0,
    name = None if channel_names is None else channel_names,
    visible = False)
    
    
    # View image
    #with napari.gui_qt():
          
    # Phenotypes under investigation
    if overlay_phenotype is None:
        available_phenotypes = list(adata.obs[phenotype_column].unique())
    else:
        available_phenotypes = overlay_phenotype
    
    # Add phenotype layer function
    def add_phenotype_layer (adata, phenotype_column, phenotype_layer,x,y,viewer,point_size,point_color):
        coordinates = adata[adata.obs[phenotype_column] == phenotype_layer]
        coordinates = pd.DataFrame({'y': coordinates.obs[y],'x': coordinates.obs[x]})
        points = coordinates.values.tolist()
        if point_color is None:
            r = lambda: random.randint(0,255) # random color generator
            point_color = '#%02X%02X%02X' % (r(),r(),r()) # random color generator
        viewer.add_points(points, size=point_size,face_color=point_color,visible=False,name=phenotype_layer)
    
    # Run the function on all phenotypes    
    for i in available_phenotypes:
        add_phenotype_layer (adata=adata, phenotype_column=phenotype_column, 
                                                           phenotype_layer=i, x=x, y=y, viewer=viewer,
                                                           point_size=point_size,point_color=point_color)
        
