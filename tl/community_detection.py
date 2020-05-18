#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 08:52:18 2020
@author: Ajit Johnson Nirmal
Nearest neighbourgh based clustering analysis

"""

from sklearn.neighbors import BallTree
import pandas as pd
import numpy as np

def community_detection (adata,x_coordinate='X_position',y_coordinate='Y_position',phenotype='phenotype',
                         neighbour_radius=30, k=5, unique_id='ImageId', nn_stat='count', method='local_radius',
                         subset_image=None,label='nn_communities'):
    """
    

    Parameters
    ----------
    adata : Ann Data Object
    x_coordinate : string, optional
        Column name of the x-coordinates. The default is 'X_position'.
    y_coordinate : string, optional
        Column name of the y-coordinates. The default is 'Y_position'.
    phenotype : string, optional
        Column name of the column containing the phenotype information. The default is 'phenotype'.
    neighbour_radius : int, optional
        The radius around each cell that is considered as a neighbhour. The default is 30.
    k : int, optional
        Number of nearest neighbours to identify. Works only when nn_stat is set to 'knn'. The default is 5.
    unique_id : string, optional
        For datasets containing multiple images, column name of the column containing the image id 
        need to be passed here. The default is 'ImageId'.
    nn_stat : string, optional
        Two options are available: a) 'count', b) 'dist'. 
        a) count- Calculates the frequency of cell phenotypes withing the cell neighbourhood
        b) dit- Calculates the mean distance of cell phrnotypes within the cell neighbourhood
        The default is 'count'.
    method : string, optional
        Two options are available: a) 'local_radius', b) 'knn'.
        a) local_radius - Indentifies the neighbours within a given radius for every cell.
        b) knn - Identifies the K nearest neigbours for every cell.
        The default is 'local_radius'.
    subset_image : string, optional
        ImageId of a single image to be subsetted for analyis. The default is None.
    label : string, optional
        Name the resulting dataframe that is returned. The default is 'nn_communities'.

    Returns
    -------
    Modified Ann Data
        Check adata.uns['nn_communities'] for results.
    
    Example
    -------
    adata = nn_communities (adata,phenotype='phenotype',neighbour_radius=30, 
                            k=5, nn_stat='count', method='local_radius',subset_image=None)


    """
    
    
    def community (adata_subset,x_coordinate,y_coordinate,phenotype,unique_id,
                   neighbour_radius, k, nn_stat, method):
    
        # Print
        print ("Analysing Image: "+ str(adata_subset.obs[unique_id].unique()))
        # Create a dataFrame with the necessary inforamtion
        data = pd.DataFrame({'x': adata_subset.obs[x_coordinate], 'y': adata_subset.obs[y_coordinate], 'phenotype': adata_subset.obs[phenotype]})
        
        # Identify neighbours
        if method == 'local_radius':
            print("Identifying neighbours within " + str(neighbour_radius) + " pixels of every cell")
            kdt = BallTree(data[['x','y']], metric='euclidean') 
            ind, dist = kdt.query_radius(data[['x','y']], r=neighbour_radius,return_distance=True)
            for i in range(0, len(ind)): ind[i] = np.delete(ind[i], np.argwhere(ind[i] == i))#remove self
            neighbours = pd.DataFrame(ind.tolist(), index = data.index) # neighbour DF
            if nn_stat == 'dist':
                for i in range(0, len(ind)): dist[i] =  dist[i][dist[i] != 0]#remove self distance
                neighbours_dist = pd.DataFrame(dist.tolist(), index = data.index) # neighbour DF          
            
        if method == 'knn':
            print("Identifying the " + str(k) + " nearest neighbours for every cell")
            tree = BallTree(data[['x','y']], leaf_size= 2)
            dist, ind = tree.query(data[['x','y']], k=k)
            neighbours = pd.DataFrame(ind.tolist(), index = data.index) # neighbour DF
            neighbours.drop(0, axis=1, inplace=True) # Remove self neighbour
            if nn_stat == 'dist':
                neighbours_dist = pd.DataFrame(dist.tolist(), index = data.index) # neighbour DF
                neighbours_dist.drop(0, axis=1, inplace=True) # Remove self neighbour
              
        # Map phenotype
        phenomap = dict(zip(list(range(len(ind))), data['phenotype'])) # Used for mapping
        
        # Loop through (all functionized methods were very slow)
        for i in neighbours.columns:
            neighbours[i] = neighbours[i].dropna().map(phenomap, na_action='ignore')
        
        # Drop NA
        #n_dropped = neighbours.dropna(how='all')
           
        # Collapse all the neighbours into a single column
        n = pd.DataFrame(neighbours.stack(), columns = ["neighbour_phenotype"])
        n.index = n.index.get_level_values(0) # Drop the multi index
        n = pd.DataFrame(n)
        n['order'] = list(range(len(n)))
        
        # Merge with real phenotype
        n_m = n.merge(data['phenotype'], how='inner', left_index=True, right_index=True)
        n_m['neighbourhood'] = n_m.index
        n = n_m.sort_values(by=['order'])
        
        if nn_stat == 'count':
            # Normalize based on total cell count
            k = n.groupby(['neighbourhood','neighbour_phenotype']).size().unstack().fillna(0)
            kraw = k
            k = k.div(k.sum(axis=1), axis=0)
        
        # Add distancce to neighbours
        if nn_stat == 'dist':
            #neighbours_dist = neighbours_dist.loc[ n_dropped.index , : ]# Drop the same in the dist DF as well
            n_dist = pd.DataFrame(neighbours_dist.stack(), columns = ["dist"])
            n_dist.index = n_dist.index.get_level_values(0) # Drop the multi index
            n_dist = pd.DataFrame(n_dist)
            n = pd.concat([n, n_dist], axis=1, sort=False)          
            # Normalize based on total cell count
            k = n.groupby(['neighbourhood','neighbour_phenotype'])['dist'].mean().unstack().fillna(0)  
            kraw = k
            k = k.div(k.sum(axis=1), axis=0)
        
        return [k, kraw]
        
        
    # Subset a particular image if needed
    if subset_image is not None:
        adata_list = [adata[adata.obs[unique_id] == subset_image]]
    else:
        adata_list = [adata[adata.obs[unique_id] == i] for i in adata.obs[unique_id].unique()]
    
    # Apply function to all images and create a master dataframe
    r_community = lambda x: community(adata_subset=x, x_coordinate=x_coordinate, y_coordinate=y_coordinate,
                                      phenotype=phenotype, unique_id=unique_id, neighbour_radius=neighbour_radius,
                                      k=k, nn_stat=nn_stat, method=method) # Create lamda function 
    all_data = list(map(r_community, adata_list)) # Apply function 


    # Merge all the results into a single dataframe    
    result = []
    for i in range(len(all_data)):
        result.append(all_data[i][0])
    result = pd.concat(result, join='outer')  
    
    result_raw = []
    for i in range(len(all_data)):
        result_raw.append(all_data[i][1])
    result_raw = pd.concat(result_raw, join='outer') 
    
    
    # Reindex the cells
    result = result.fillna(0)
    result = result.reindex(adata.obs.index)
    
    result_raw = result_raw.fillna(0)
    result_raw = result_raw.reindex(adata.obs.index)
    
    # Add to adata
    adata.uns[label] = result
    adata.uns[str(label)+"_raw"] = result_raw
    
    # Also add to adata.obs
    for i in result.columns:
        adata.obs[i] = result[i]
        
    # Return
    return adata

    