# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 13:54:38 2020
@author: Ajit Johnson Nirmal
Using the mutually exclusive markers, perform background normalization (Implementation of RESTORE method)
"""
# Import library
import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture 
from sklearn.cluster import KMeans
import seaborn as sns; sns.set(color_codes=True)
import matplotlib.pyplot as plt
from preprocess.selfrepresentation import SparseSubspaceClusteringOMP


# Background removal
def background_normalization (adata, mutually_exclusive_markers, replace_marker=None, seed=100):
    '''
    Parameters:
        adata- AnnData object created using the function df_to_annobject 
        mutually_exclusive_markers: Pandas dataframe containing the mutually exclusive markers; either manually curated or determined using mutually_exclusive_markers function.
        seed: int. Random seed.
        replace_marker: Nested list. Use this to replace any mutually exclusive marker pairs.
    Returns:
        AnnData with normlized data
    Example:
        background_normalization(adata, mutually_exclusive_markers, replace_marker=[['CD11B', 'ASMA'],['CD11C', 'ASMA']], seed=100)
    '''
     
    # Prepare the mutually exclusive dataframe (subset the first two columns)
    exclusive_markers = mutually_exclusive_markers.iloc[:,:2]
    exclusive_markers.columns = ['MA', 'MB']
    if replace_marker != None:
        # Replace manual elements
        for i in range(len(replace_marker)):
            if  any(exclusive_markers['MA'].str.contains(replace_marker[i][0])):
                exclusive_markers['MB'][np.where(exclusive_markers['MA'].str.contains(replace_marker[i][0]))[0]] = replace_marker[i][1]
    # Modified marker list to list
    exclusive_markers = exclusive_markers.values.tolist()
 
    # Function to cluster
    def clustering (data, marker_pair):
        
        # Subset data
        data = pd.DataFrame(list(zip(adata[:,marker_pair[0]].X, adata[:,marker_pair[1]].X)), columns =[marker_pair[0], marker_pair[1]])

        # Cluatering
        scc = SparseSubspaceClusteringOMP(n_clusters = 2, random_state = seed).fit(data)
        gmm = GaussianMixture(n_components = 2, random_state=seed).fit(data).predict(data)
        km = KMeans(n_clusters = 2, random_state=seed).fit(data).predict(data)
        
        # Add the clusters as a column to the data
        data['scc'] = scc.labels_
        data['gmm'] = gmm
        data['km'] = km
        
        # Find the negative cluster and get the mean and STD of the negative cluster
        # SCC
        if np.mean(data[data['scc'] == 0].iloc[:,0]) < np.mean(data[data['scc'] == 1].iloc[:,0]):
            scc_c = [np.mean(data[data['scc'] == 0].iloc[:,0]), np.std(data[data['scc'] == 0].iloc[:,0])]
        else:
            scc_c = [np.mean(data[data['scc'] == 1].iloc[:,0]), np.std(data[data['scc'] == 1].iloc[:,0])]
        
        # GMM
        if np.mean(data[data['gmm'] == 0].iloc[:,0]) < np.mean(data[data['gmm'] == 1].iloc[:,0]):
            gmm_c = [np.mean(data[data['gmm'] == 0].iloc[:,0]), np.std(data[data['gmm'] == 0].iloc[:,0])]
        else:
            gmm_c = [np.mean(data[data['gmm'] == 1].iloc[:,0]), np.std(data[data['gmm'] == 1].iloc[:,0])]
            
        #KM
        if np.mean(data[data['km'] == 0].iloc[:,0]) < np.mean(data[data['km'] == 1].iloc[:,0]):
            km_c = [np.mean(data[data['km'] == 0].iloc[:,0]), np.std(data[data['km'] == 0].iloc[:,0])]
        else:
            km_c = [np.mean(data[data['km'] == 1].iloc[:,0]), np.std(data[data['km'] == 1].iloc[:,0])]
        
        # Calculate the median of the three clustering methods
        mean = np.mean([scc_c[0], gmm_c[0], km_c[0]])
        std = np.mean([scc_c[1], gmm_c[1], km_c[1]])
        
        # Background intensity
        bi = mean + std
        
        # Substract the background intensity from the ith marker
        data_normalized = data.iloc[:,0] - bi
        data_normalized[data_normalized < 0] = 0 # Remove less than zero values
               
        # Return data
        return data_normalized
       
    # Run the clustering function on all markers (vectorized)
    r_clustering = lambda x: clustering(marker_pair = x, data = adata) # Create lamda function 
    normalized_data = list(map(r_clustering, exclusive_markers)) # Apply function
    normalized_data = pd.concat(normalized_data, axis=1, keys=[s.name for s in normalized_data])
    
    # create a copy of unnormalized data
    adata.raw = adata
    
    # Replace with normalized data
    normalized_data = normalized_data[adata.var.index]
    adata.X = np.array(normalized_data)
    
    # retun data
    return adata
    
