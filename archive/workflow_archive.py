# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 13:37:55 2020
@author: ajn16
WorkFlow Archive
"""

from preprocess.mutually_exclusive_markers import mutually_exclusive_markers
from preprocess.selfrepresentation import SparseSubspaceClusteringOMP
from preprocess.background_modeling import background_modeling
from plots.raw_vs_normalized_plot import raw_vs_normalized_plot

# Identify the mutually exclusive markers for background removal
adata = mutually_exclusive_markers (adata, top_features=3, log=True, random_sample=5000)
mutually_exclusive_markers = pd.read_csv("C:/Users/ajn16/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/3- background removal/PTCL1/mutually_exclusive_markers.csv")

# View results
mutually_exclusive_markers = adata.uns['mutually_exclusive_markers']
#mutually_exclusive_markers.to_csv("C:/Users/ajn16/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/3- background removal/PTCL1/mutually_exclusive_markers.csv")
adata.uns['mutually_exclusive_markers_with_scores']
svd_ratio= adata.uns['svd_ratio_matrix'] 
# Heatmap of the SVD
sns.clustermap(svd_ratio)


# Remove background to identify positivity for each marker

replace_marker = [['CD11B', 'ASMA'],['CD11C', 'ASMA'],
                  ['CD20', 'ASMA'],['FOXP3', 'ASMA'],
                  ['HLADR', 'ASMA'],['KI67', 'ASMA'],
                  ['PS6', 'ASMA']]

out_dir= 'C:/Users/ajn16/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/3- background removal/PTCL1/prediction_plots'
adata = background_modeling(adata, mutually_exclusive_markers, sample_data=True, 
                            fraction=0.1, random_sample=None, seed=100, 
                            replace_marker=None, plot_figure=True, out_dir=out_dir, 
                            figure_type="png", percentile=[75,99])


# Generate plots for comparision of before and after normalization
n_plots = raw_vs_normalized_plot (adata, out_dir=out_dir,log=True)


