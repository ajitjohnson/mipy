# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 12:13:22 2020
@author: Ajit Johnson Nirmal
ReScale the data
"""
# Import library
import pandas as pd
import numpy as np
import itertools
from sklearn.preprocessing import MinMaxScaler
from sklearn.mixture import GaussianMixture 


def rescale (adata, gate=None, return_gates=False, failed_markers=None, method='all'):
    '''
    Parameters:
        data: AnnData object
        gate: DataFrame with first column as markers and second column as the gate values in log1p scale
        failed_markers: list. list of markers that are not expressed at all in any cell. pass in as ['CD20', 'CD3D']
    Returns:
        Ann Data object with rescaled data
    Example:
        adata = rescale (adata, gate=manual_gate, failed_markers=['CD20', 'CD21'])
    '''
    def rescale_independent (adata, gate=None, return_gates=False, failed_markers=None):
        
        print('Scaling Image '+ str(adata.obs['ImageId'].unique()))
        
        # Copy of the raw data if it exisits
        if adata.raw is not None:
            adata.X = adata.raw.X
              
        data = pd.DataFrame(adata.X, columns = adata.var.index, index= adata.obs.index)      
        # Merging the manual gates and non-working markers togeather if any
        if gate is not None:
            m_markers = list(gate.iloc[:,0])
            manual_gate_markers = gate
        if failed_markers != None:
            manual_gate_markers = pd.DataFrame(data[failed_markers].quantile(0.9999999))
            manual_gate_markers['markers'] = failed_markers
            # move column to front
            cols = manual_gate_markers.columns.tolist()
            cols.insert(0, cols.pop(cols.index('markers')))
            manual_gate_markers = manual_gate_markers.reindex(columns= cols)
            manual_gate_markers.columns = ['marker', 'gate']
            m_markers = failed_markers
        if gate is not None and failed_markers != None:
            m_markers = list(gate.iloc[:,0]) + list(manual_gate_markers.iloc[:,0])
            gate.columns = ['marker', 'gate']
            manual_gate_markers = pd.concat([gate, manual_gate_markers])
        if gate is None and failed_markers == None:
            m_markers = []
                  
        # Find markers to send to gmm modelling
        if gate is not None or failed_markers is not None:
            gmm_markers = list(np.setdiff1d(data.columns, m_markers))
        else:
            gmm_markers = list(data.columns)
            
        # If manual gate is not provided scale the data 
        if len(gmm_markers) != 0:
            gmm_data = data[gmm_markers]        
            # Clip off the 99th percentile    
            def clipping (x):
                clip = x.clip(lower =np.percentile(x,1), upper=np.percentile(x,99)).tolist()
                return clip
            # Run the function
            gmm_data = gmm_data.apply(clipping)
            # Scaling the data
            sum_data = gmm_data.sum(axis=1) # Calculate total count for each cell
            n_count = gmm_data.div(sum_data, axis=0) # Divide genes by total count for every cell        
            med = np.median(list(itertools.chain(*gmm_data.values.tolist()))) # Calculate median count of the entire dataset
            n_count = n_count*med # Multiply by scaling fator (median count of entire dataset)
            n_log = np.log1p(n_count) # Log transform data
            scaler = MinMaxScaler(feature_range=(0, 1))
            s = scaler.fit_transform(n_log)
            normalised_data = pd.DataFrame(s, columns = gmm_data.columns, index= gmm_data.index)
            # Gaussian fit to identify the gate for each marker and scale based on the gate
            # Empty data frame to hold the results
            all_gmm_data = pd.DataFrame()
            def gmm_gating (data, marker, return_gates):
                # Print
                print('Finding the optimal gate for ' + str(marker))
                # Identify the marker to fit the model
                m = data[marker].values
                # Perform GMM
                data_gm = m.reshape(-1, 1)
                #gmm = GaussianMixture(n_components=2, means_init=[[0],[1]],covariance_type='tied')
                gmm = GaussianMixture(n_components=2)
                gmm.fit(data_gm)
                gate = np.mean(gmm.means_)
                
                # Find the closest value to the gate
                absolute_val_array = np.abs(m - gate)
                smallest_difference_index = absolute_val_array.argmin()
                closest_element = m[smallest_difference_index]
                
                # rescale the data based on the identified gate
                marker_study = pd.DataFrame(m, index= data.index)
                marker_study = marker_study.sort_values(0)
                
                # Find the index of the gate
                gate_index = marker_study.index[marker_study[0] == closest_element][0]
                
                # Split into high and low groups 
                high = marker_study.loc[gate_index:,:]
                low = marker_study.loc[:gate_index,:]
                
                # Prepare for scaling the high and low dataframes
                scaler_high = MinMaxScaler(feature_range=(0.5, 1))
                scaler_low = MinMaxScaler(feature_range=(0, 0.5))
                
                # Scale it
                h = pd.DataFrame(scaler_high.fit_transform(high), index = high.index)
                l = pd.DataFrame(scaler_low.fit_transform(low), index = low.index)
                
                # Merge the high and low and resort it
                scaled_data = pd.concat([l,h])
                scaled_data = scaled_data.loc[~scaled_data.index.duplicated(keep='first')]           
                scaled_data = scaled_data.reindex(data.index)
    
                #return scaled_data
                if return_gates == True:
                    return gate
                else:
                    return scaled_data
            
            # Apply the function
            r_gmm_gating = lambda x: gmm_gating(data=normalised_data, marker=x,return_gates=return_gates) # Create lamda function 
            all_gmm_data = list(map(r_gmm_gating, gmm_markers)) # Apply function        
            all_gmm_data = pd.concat(all_gmm_data, axis=1, sort=False)        
            all_gmm_data.columns = gmm_markers
        else:
            all_gmm_data = pd.DataFrame()
        
        # Empty data frame to hold the results
        all_manual_data = pd.DataFrame()
        if len(m_markers) != 0:
            m_data = np.log1p(data[m_markers]) 
            # Clip the data
            def clipping (x):
                clip = x.clip(lower =np.percentile(x,1), upper=np.percentile(x,99)).tolist()
                return clip
            # Run the function
            m_data = m_data.apply(clipping)
            
            def manual_gating (data,marker,gate):
                # Print
                print('Scaling ' + str(marker))
                # Work on processing manual gates
                m = gate[gate.iloc[:,0] == marker].iloc[:,1] # gate of the marker passed in
                
                # Find the closest value to the gate
                absolute_val_array = np.abs(data[marker].values - float(m))
                smallest_difference_index = absolute_val_array.argmin()
                closest_element = data[marker].values[smallest_difference_index]
                
                # rescale the data based on the identified gate
                marker_study = data[marker]
                marker_study = marker_study.sort_values(0)
                
                # Find the index of the gate
                gate_index = marker_study.index[marker_study == closest_element][0]
                
                # Split into high and low groups 
                high = marker_study[gate_index:]
                low = marker_study[:gate_index]
                
                # Prepare for scaling the high and low dataframes
                scaler_high = MinMaxScaler(feature_range=(0.5, 1))
                scaler_low = MinMaxScaler(feature_range=(0, 0.5))
                
                # Scale it
                h = pd.DataFrame(scaler_high.fit_transform(high.values.reshape(-1, 1)), index = high.index)
                l = pd.DataFrame(scaler_low.fit_transform(low.values.reshape(-1, 1)), index = low.index)
                
                # Merge the high and low and resort it
                scaled_data = pd.concat([l,h])
                scaled_data = scaled_data.loc[~scaled_data.index.duplicated(keep='first')]           
                scaled_data = scaled_data.reindex(data.index)
                
                # Return
                return scaled_data
            
            # Apply the function
            r_manual_gating = lambda x: manual_gating(data=m_data, marker=x, gate=manual_gate_markers) # Create lamda function 
            all_manual_data = list(map(r_manual_gating, m_markers)) # Apply function        
            all_manual_data = pd.concat(all_manual_data, axis=1, sort=False)        
            all_manual_data.columns = m_markers
        
        else:
            all_manual_data = pd.DataFrame()
            
            
        # If both manual and automatic gating was used, combine them into a single result    
        if not all_manual_data.empty:
            all_scaled_data = all_manual_data
        if not all_gmm_data.empty:
            all_scaled_data = all_gmm_data
        if not all_manual_data.empty and not all_gmm_data.empty:
            all_scaled_data = all_gmm_data.merge(all_manual_data, how='outer', left_index=True, right_index=True)
        
        # re index the columns
        all_scaled_data = all_scaled_data.reindex(columns= data.columns)
        return all_scaled_data
    
    # Apply method of choice
    if method == 'all':
        all_scaled_data = rescale_independent (adata, gate=gate, return_gates=return_gates, failed_markers=failed_markers)
    if method == 'by_image':
        adata_list = [adata[adata.obs['ImageId'] == i] for i in adata.obs['ImageId'].unique()]
        r_rescale_independent = lambda x: rescale_independent(adata=x, gate=gate, return_gates=return_gates, failed_markers=failed_markers) # Create lamda function 
        scaled_data = list(map(r_rescale_independent, adata_list)) # Apply function
        all_scaled_data = pd.concat(scaled_data)
                   
    # Create a copy of the raw data
    if adata.raw is None:
        adata.raw = adata
     
    # Replace with normalized data
    adata.X = all_scaled_data
    
    # Return data
    return adata

