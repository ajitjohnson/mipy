# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 13:54:38 2020
@author: Ajit Johnson Nirmal
Using the mutually exclusive markers, perform background normalization (Implementation of modified RESTORE method)
"""
# Library
import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture 
from sklearn.cluster import KMeans
import seaborn as sns; sns.set(color_codes=True)
import matplotlib.pyplot as plt
import scanpy as sc
# Modelling
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import accuracy_score
from sklearn import metrics

def background_modeling (adata, mutually_exclusive_markers, percentile=[75,99], sample_data=True, fraction=None, random_sample=None, seed=100, replace_marker=None, plot_figure=True, out_dir=None, figure_type="png"):
    '''
    Parameters:
        adata- AnnData object created using the function df_to_annobject 
        mutually_exclusive_markers: Pandas dataframe containing the mutually exclusive markers; either manually curated or determined using mutually_exclusive_markers function.
        sample_data: Boolian. If true, either use 'fraction' or 'random_sample' to randomly choose a subset of the dataset.
        percentile: List of two int. The lower and upper bound for choosing positive cells.
        fraction: use valuses between [0-1]. Subsample to this fraction of the number of observations.
        random_sample: int. Subsample to this number of observations.
        seed: int. Random seed.
        replace_marker: Nested list. Use this to replace any mutually exclusive marker pairs.
        plot_figure: Boolian. If True, scatter plot for every marker showing cells that were chosen as positive control will be displayed.
        out_dir: Directory to which the figures will be saved. If no directory is passed, figures will be saved in the working directory.
        figure_type: Png or Pdf.
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
    
    # Randomly sample data
    if sample_data == True:
        bdata = sc.pp.subsample(adata,fraction=fraction, n_obs=random_sample, random_state=seed, copy=True)
    else:
        bdata = adata
    
    # Full dataset to pass in for prediction as a dataframe
    full_data = pd.DataFrame(adata.X, columns = adata.var.index, index= adata.obs.index)
    #full_data = np.log1p(full_data)
        
    
    def predict_positivity (data, full_data, marker_pair, plot_figure, out_dir, figure_type, percentile):
        print(marker_pair)
        # Subset data
        data = pd.DataFrame(list(zip(data[:,marker_pair[0]].X, data[:,marker_pair[1]].X)), columns =[marker_pair[0], marker_pair[1]])
        data.index = bdata.obs.index
        # Clustering
        gmm = GaussianMixture(n_components = 2, random_state=seed).fit(data).predict(data)
        km = KMeans(n_clusters = 2, random_state=seed).fit(data).predict(data)
        # Add the clusters as a column to the data
        data['gmm'] = gmm
        data['km'] = km
        ## NEGATIVE CLUSTER- mean and STD of i-th marker
        # GMM
        if np.mean(data[data['gmm'] == 0].iloc[:,0]) < np.mean(data[data['gmm'] == 1].iloc[:,0]):
            gmm_n = [np.mean(data[data['gmm'] == 0].iloc[:,0]), np.std(data[data['gmm'] == 0].iloc[:,0]), np.percentile(data[data['gmm'] == 0].iloc[:,1] , 25)]
        else:
            gmm_n = [np.mean(data[data['gmm'] == 1].iloc[:,0]), np.std(data[data['gmm'] == 1].iloc[:,0]), np.percentile(data[data['gmm'] == 1].iloc[:,1] , 25)]
            
        #KM
        if np.mean(data[data['km'] == 0].iloc[:,0]) < np.mean(data[data['km'] == 1].iloc[:,0]):
            km_n = [np.mean(data[data['km'] == 0].iloc[:,0]), np.std(data[data['km'] == 0].iloc[:,0]), np.percentile(data[data['km'] == 0].iloc[:,1] , 25)]
        else:
            km_n = [np.mean(data[data['km'] == 1].iloc[:,0]), np.std(data[data['km'] == 1].iloc[:,0]), np.percentile(data[data['km'] == 1].iloc[:,1] , 25)]
        
        ## POSITIVE CLUSTER- 75th and 99th percentile percentile[1]
        # IKM
        if np.mean(data[data['km'] == 0].iloc[:,0]) < np.mean(data[data['km'] == 1].iloc[:,0]):
            km_p = [np.percentile(data[data['km'] == 1].iloc[:,0] , percentile[0]), np.percentile(data[data['km'] == 1].iloc[:,0], percentile[1])]
        else:
            km_p = [np.percentile(data[data['km'] == 0].iloc[:,0], percentile[0]), np.percentile(data[data['km'] == 0].iloc[:,0], percentile[1])]
            
        # GMM
        if np.mean(data[data['gmm'] == 0].iloc[:,0]) < np.mean(data[data['gmm'] == 1].iloc[:,0]):
            gmm_p = [np.percentile(data[data['gmm'] == 1].iloc[:,0] , percentile[0]), np.percentile(data[data['gmm'] == 1].iloc[:,0], percentile[1])]
        else:
            gmm_p = [np.percentile(data[data['gmm'] == 0].iloc[:,0], percentile[0]), np.percentile(data[data['gmm'] == 0].iloc[:,0], percentile[1])]   
        
        # Mean of the Negative cluster
        # Calculate the median of the three clustering methods
        negative_mean = np.mean([gmm_n[0], km_n[0]])
        negative_std = np.mean([gmm_n[1], km_n[1]])
        negative = negative_mean + negative_std
        # Upper threshold/ negative leak into postive cells
        upper_threshold = np.mean([gmm_n[2], km_n[2]])    
        # Calculate the mean 75th and 95th percentile from both clustering methods
        positive_75 = np.mean([gmm_p[0], km_p[0]])
        positive_99 = np.mean([gmm_p[1], km_p[1]])
        positive = np.max([negative, positive_75]) # Identify the starting point of positivity
        # Based on the defined 75th and 99th percentile, identify the positive cells
        positive_cells = data[data.iloc[:,0].between(positive, positive_99)]
        positive_cells = positive_cells[positive_cells.iloc[:,1] < upper_threshold]
        
        ############################
        # Model the data
        ############################
        # Get the real data
        real_data = pd.DataFrame(bdata.X, columns = bdata.var.index, index = bdata.obs.index)
        #real_data = np.log1p(real_data)
        # Add the positive cells label to the df
        real_data['label'] = np.where(real_data.index.isin(positive_cells.index), 1, 0)
        # Create a new df with equal representation of positive and negative cells
        no_of_cells = np.min([len(real_data[real_data['label'] == 1]), len(real_data[real_data['label'] == 0])])
        # New dataframe to be used for modelling
        data_model = real_data[real_data['label'] == 0].sample(n=no_of_cells).append(real_data[real_data['label'] == 1].sample(n=no_of_cells))
        data_model = data_model.sample(frac=1) # Shuffle the dataframe
        
        # Create the train-test split
        train, test = train_test_split(data_model, test_size = 0.4)
        train = train.reset_index(drop=True)
        test = test.reset_index(drop=True)
        
        # Features
        features_train = train[train.columns.drop('label')]
        label_train = train['label']
        features_test = test[test.columns.drop('label')]
        label_test = test['label']
        
        # Preprocess the training and testing data
        scaler = StandardScaler()
        scaler.fit(features_train)
        features_train = scaler.transform(features_train)
        
        # Apply same trasformation to test data
        features_test = scaler.transform(features_test)  
        
        # Build the Model
        clf = MLPClassifier()
        clf.fit(features_train,label_train)
        pred_train = clf.predict(features_train)
        pred_test = clf.predict(features_test)
        accuracy_train = accuracy_score(pred_train,label_train)
        accuracy_test = accuracy_score(pred_test,label_test)
        fpr, tpr, _ = metrics.roc_curve(np.array(label_train), clf.predict_proba(features_train)[:,1])
        auc_train = metrics.auc(fpr,tpr)
        fpr, tpr, _ = metrics.roc_curve(np.array(label_test), clf.predict_proba(features_test)[:,1])
        auc_test = metrics.auc(fpr,tpr)
        #print(accuracy_train,accuracy_test,auc_train,auc_test)
        
        # Prediction on the whole dataset
        f = scaler.transform(full_data)  
        prediction = clf.predict(f)
        
        # Figure
        if plot_figure == True:
            prediction_accuracy = "Prediction Accuracy: " + str(round(accuracy_test, 2))
            prediction_auc = "Prediction AUC: " + str(round(auc_test, 2))
            No_cells = "No. of cells in +ve bin: " + str(len(positive_cells))
            # Figure
            sns.set_style("white")
            sns.scatterplot(x=data.iloc[:,0], y=data.iloc[:,1],hue=km) 
            plt.title(str(marker_pair) + ' - K Means Clustering')
            plt.plot([positive, positive], [0, upper_threshold], color= 'k')
            plt.plot([positive_99, positive_99], [0, upper_threshold], color= 'k')
            plt.plot([positive, positive_99], [upper_threshold, upper_threshold], color= 'k')
            plt.plot([positive, positive_99], [0, 0], color= 'k')
            plt.annotate(No_cells, (0,0), (0, -35), xycoords='axes fraction', textcoords='offset points', va='top')
            plt.annotate(prediction_accuracy, (0,0), (0, -45), xycoords='axes fraction', textcoords='offset points', va='top')
            plt.annotate(prediction_auc, (0,0), (0, -55), xycoords='axes fraction', textcoords='offset points', va='top')
            if out_dir != None:
                plt.savefig(out_dir + "/" + marker_pair[0] + "." + figure_type, bbox_inches="tight")
            else:
                plt.savefig(marker_pair[0] + "." + figure_type, bbox_inches="tight")  
            plt.clf()
            
        # Return the predictions
        return prediction
    

    # Run the clustering function on all markers (vectorized)
    r_predict = lambda x: predict_positivity(marker_pair = x, data = bdata, full_data=full_data, plot_figure=plot_figure, out_dir=out_dir, figure_type=figure_type, percentile=percentile) # Create lamda function 
    predictions = list(map(r_predict, exclusive_markers)) # Apply function
    predictions = pd.DataFrame(np.row_stack(predictions)).T.astype('category')
    predictions.columns = list(next(zip(*exclusive_markers)))
    predictions.index = adata.obs.index
    
    # Return data
    obs = adata.obs
    result = pd.concat([obs, predictions], axis=1, sort=False)
    adata.obs = result
    
    # Also save in the uns region
    adata.uns['prediction'] = predictions
    
    # return adata
    return adata
