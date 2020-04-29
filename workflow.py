# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 09:12:39 2020
@author: Ajit Johnson Nirmal
MiPy package functions workflow- PTCL TMA
"""

# Import necesaary packages
import sys
import os
import pandas as pd
import scanpy as sc
import seaborn as sns; sns.set(color_codes=True)
import scanpy.external as sce
import matplotlib.pyplot as plt

# Import data
adata = sc.read('/Users/aj/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/4-celltype_calling/PTCL6/PTCL6.h5ad')


# Set WD
os.chdir("/Users/aj/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/3- background removal")

#  Import SiPy package components
sys.path.insert(1, '/Users/aj/Dropbox (Partners HealthCare)/packages/mipy')
sys.path.insert(1, 'C:/Users/ajn16/Dropbox (Partners HealthCare)/packages/mipy')
from preprocess.mcmicro_to_object import mcmicro_to_object
from plots.marker_expression_plot import marker_expression_plot
from preprocess.rescale import rescale
from phenotyping.phenotype_groups import phenotype_groups
from phenotyping.phenotype_subcluster import phenotype_subcluster
from processing.rename_clusters import rename_clusters
from spatial.knn import knn
from spatial.local_neighbourhood import local_neighbourhood
from plots.napari import *
from export.save_files import *
#===============================================================================

# Import data
#data = pd.read_csv('/Users/aj/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/celltype_calling/PTCL1.csv')
#data = pd.read_csv('/Users/ajn16/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/celltype_calling/PTCL1_subset.csv')
#data = pd.read_csv('//research.files.med.harvard.edu/ImStor/sorger/data/RareCyte/Connor/Z155_PTCL/whole_section_csv_files/raw/PTCL1_450.csv')

# Preporcess the data
adata = histocat_to_object (data, remove_string_from_name = "Cell_PTCL1_450", islog=True, drop_markers= ["PERK", "NOS2"])

# Processing multiple images
image_path = ['/Users/aj/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/raw_data/whole_sections/PTCL1_450.csv',
              '/Users/aj/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/raw_data/whole_sections/PTCL2_552.csv']

adata = histocat_to_object (image_path, remove_string_from_name = "Cell_PTCL1_450", islog=True, drop_markers= ["PERK", "NOS2"])

sc.pp.combat(adata, key='ImageId')


import numpy as np
adata.X = np.exp(adata.X)


# To reuse adata
for i in adata.var.index:
    try:
        del adata.obs[i]      
    except:
        pass
# To reuse adata
        
for i in phenotype.iloc[:,0].unique():
    try:
        del adata.obs[i]      
    except:
        pass
    
# Rescale data
manual_gates = pd.read_csv('/Users/aj/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/4-celltype_calling/markers copy.csv')

adata.X = adata.raw.X
# Without manual gates
#adata = rescale (adata, gate=None, failed_markers=['CD20','CD21','CD56','CD15'])
adata = rescale (adata, gate=None)
# With manual gates
adata = rescale (adata, gate=manual_gates, return_gates=False)

# Create a positivity for all markers
out_dir= '/Users/aj/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/4-celltype_calling/PTCL2/marker_expression_plot'
plots = marker_expression_plot (adata, x='X_position', y='Y_position', 
                                    color=None, save_figure= True, plot_all=True, 
                                    out_dir=out_dir, alpha= 0.3, figure_type="png", 
                                    s=1, cmap='gist_heat', figsize=(12, 5))

# Phenotyping
phenotype = pd.read_csv('/Users/aj/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/4-celltype_calling/PTCL1/phenotype_workflow.csv')
adata = phenotype_groups (adata, phenotype, gate = 0.5, label="phenotype")
# summary stats
adata.obs['phenotype'].value_counts()

# Rename clusters
rename = {'other non-immune cells': ['non-immune'],
          'other immune cells': ['immune']}
adata = rename_clusters (adata, rename, from_column='phenotype', to_column='phenotype_renamed')
# summary stats
adata.obs['phenotype_renamed'].value_counts()

# Sub-phenotyping
adata = phenotype_subcluster (adata, cluster_group = ['tumor'], k = 15, method = 'kmeans',label='phenotype_renamed', use_raw_data = True)
adata.obs['kmeans'].value_counts()

# Tumor marker genes
tumor = adata[adata.obs['phenotype_renamed'] == 'tumor']
sc.tl.rank_genes_groups(tumor, groupby= 'kmeans',n_genes=5)
sc.pl.rank_genes_groups(tumor, n_genes=5, sharey=False, size = 40, fontsize= 20)
sc.pl.matrixplot(tumor, var_names= tumor.var.index, groupby='kmeans', dendrogram=True, use_raw=True, cmap=c, standard_scale='var', log=True)
sc.pl.matrixplot(tumor, var_names= tumor.var.index, groupby='kmeans', dendrogram=True, use_raw=False, cmap=c,vmin=0,vmax=1,standard_scale='var')
sc.pl.matrixplot(tumor, var_names= tumor.var.index, groupby='kmeans', dendrogram=True, use_raw=False, cmap='viridis_r',vmin=0.5,vmax=1)

# Rename clusters
rename = {'tumor- KI67+': ['tumor-4', 'tumor-8', 'tumor-10'],
          'tumor- FOXP3+': ['tumor-14'],
          'tumor- PD1 high': ['tumor-1','tumor-11','tumor-12','tumor-5','tumor-2','tumor-6','tumor-13'],
          'ASMA+ cells': ['tumor-9'],
          'tumor- PD1 low': ['tumor-0','tumor-7','tumor-3']}

adata = rename_clusters (adata, rename, from_column='kmeans', to_column='kmeans_renamed')
adata.obs['kmeans_renamed'].value_counts()

# Marker genes for each cluster
sc.tl.rank_genes_groups(adata, groupby= 'kmeans_renamed',n_genes=5)
sc.pl.rank_genes_groups(adata, n_genes=5, sharey=False, size = 40, fontsize= 20)
pd.DataFrame(adata.uns['rank_genes_groups']['names'])

# PAGA- tumor
tumor = adata[adata.obs['phenotype_renamed'] == 'tumor']
sc.pp.neighbors(tumor, n_neighbors=10)
sc.tl.paga(tumor, groups='kmeans_renamed')
sc.pl.paga(tumor,layout='fr',threshold=0.2,fontsize=15,node_size_scale=5,node_size_power=0.1, max_edge_width=1)

# Non-Tumor marker genes
ntumor = adata[adata.obs['phenotype_renamed'] != 'tumor']
sc.tl.rank_genes_groups(ntumor, groupby= 'kmeans_renamed',n_genes=5)
sc.pl.rank_genes_groups(ntumor, n_genes=5, sharey=False, size = 40, fontsize= 20)
sc.pl.matrixplot(ntumor, var_names= ntumor.var.index, groupby='kmeans_renamed', dendrogram=True, use_raw=True, cmap=c, standard_scale='var', log=True)

# PAGA- nontumor
sc.pp.neighbors(ntumor, n_neighbors=10)
sc.tl.paga(ntumor, groups='kmeans_renamed')
sc.pl.paga(ntumor,layout='fa',threshold=0.4,fontsize=10,node_size_scale=3,node_size_power=0.2, max_edge_width=1)

# Save Phenotype
out_dir = '/Users/aj/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/4-celltype_calling/PTCL1/'
save_phenotype (adata,phenotype_column='kmeans_renamed',file_name='PTCL phenotype',out_dir=out_dir)
# Save Phenotype proportion
save_phenotype_proportion (adata,phenotype_column='kmeans_renamed',file_name='ptcl1_phenotype_proportion',out_dir=out_dir)

# View images
# CSV file
adata = sc.read('/Users/aj/Desktop/exemplar-002/examplar4.h5ad')
image_path = '/Users/aj/Desktop/exemplar-002/dearray/4.tif'
image_path = '/Users/aj/Desktop/exemplar-001/registration/exemplar-001.ome.tif'
I = I[0, :, 0, :, :]
view_image (image_path, adata, phenotype_column='phenotype', 
                overlay_phenotype=['Keratin','cytotoxic t cells'],
                markers=['KERATIN','CD11B', 'SMA', 'CD16'],
                channel_names=adata.uns['all_markers'],
                x='X_position',y='Y_position',point_size=10)



# UMAP of a subset of the cells
bdata = sc.pp.subsample(adata, fraction=0.1, copy=True)
sc.settings.n_jobs=4
sc.tl.pca(bdata)
sc.pp.neighbors(bdata, n_neighbors=10)
sc.tl.umap(bdata)
sc.pl.umap(bdata, color='kmeans_renamed', size = 15, use_raw=False, color_map= c, legend_loc='on data')

# HeatMap of the defined clusters
markers = ['CD45','CD3D','CD2','CD5','CD7','KI67','FOXP3','CD11C','CD11B',
           'CD163','CD206','CD68','HLADR','CD21','PDL1','PD1','CD25','CD8A','CD4','CD20','ASMA']
sc.pl.heatmap(adata, markers, groupby='kmeans_renamed', swap_axes=True, use_raw=False, log=False, cmap= 'vlag')
sc.pl.heatmap(adata, markers, groupby='kmeans_renamed', swap_axes=True, use_raw=True, log=True, cmap= 'vlag')

# Nerest Neighbour analysis
adata = knn (adata,phenotype='kmeans_renamed',permutation=1000,p_val=0.001)
neighbours_max = adata.uns['knn']

# Plot nearest neighbors
mask = neighbours_max.isnull()   
sns.set_style("whitegrid", {'axes.grid' : True})
sns.clustermap(neighbours_max.fillna(0), cmap = 'vlag', mask=mask, row_cluster=False, col_cluster=False)



# create custom color pallate
flatui = ["#515eac", "#7BCF9F", "#F7F8B7", "#FFB055", "#F64641", "#B00042"]
sns.set_palette(flatui)


import numpy as np

sns.set(rc={'figure.figsize':(11.7,8.27)})
# cluster plot not working
a = np.where(adata.obs['phenotype'] == 'Non-immune cells', 1, 0)
c_col = ['red' if x==1 else 'black' for x in a]
sns.set_style("white")
fig, ax = plt.subplots()
plt.scatter(x=adata.obs['X_position'], 
            y=adata.obs['Y_position'], 
            c=c_col, edgecolors='none', s=2, alpha=0.8)
ax.invert_yaxis()

adata.obs['phenotype'].value_counts()


geneA = 'PD1'
a = np.where(adata[:,geneA].X >= 0.5, 1, 0)
c_col = ['red' if x==1 else 'black' for x in a]
sns.set_style("white")
fig, ax = plt.subplots()
plt.scatter(x=adata.obs['X_position'], 
            y=adata.obs['Y_position'], 
            c=c_col, edgecolors='none', s=2, alpha=0.8)
ax.invert_yaxis()




# plot the size of the cells

geneA = 'CD45'
#full = adata[:,geneA].X
#full = np.log1p(adata.raw[:,geneA].X)
sns.distplot(full, rug=False, hist=True)



# Multi gate
B_T_overlap = all_scaled_data[ (all_scaled_data['CD45'] > cutoff) & (all_scaled_data['CD20'] > cutoff) & (all_scaled_data['CD3D'] > cutoff)]


# PCA variance ratio
sc.tl.pca(tumor_cells, svd_solver='arpack')
sc.pl.pca(tumor_cells, color='CD3D')
sc.pl.pca_variance_ratio(tumor_cells, log=True)

# Neighbourhood graph
sc.pp.neighbors(tumor_cells, n_neighbors=30, n_pcs=20)

# UMPA
color_map = plt.cm.get_cmap('RdBu')
c = color_map.reversed()


sc.tl.umap(tumor_cells)
sc.pl.umap(tumor_cells, color=['FOXP3', 'PSTAT3', 'PS6'], use_raw= False, size = 50, color_map= c)
sc.pl.umap(tumor_cells, color=[CD20', 'CD3D', 'CD163'], use_raw= True, size = 50, color_map= 'gist_heat')

# Palantir
d = sce.tl.palantir(adata)


#
df = pd.DataFrame(adata.var.index)

scatter_plot (adata, exp =  np.log1p(adata.X), x='X_position', y='Y_position', marker= 'CD45')
scatter_plot (adata, exp = np.log1p(adata.raw.X), x='X_position', y='Y_position', marker= 'CD45')

sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.set_style("white")
cluster_plot (adata, color ='CD45', x='X_position', y='Y_position',  s=0.1, cmaps='gist_heat')

adata

sc.pp.subsample(adata, n_obs=10000)


data = pd.read_csv(image_path[0])
data.columns

del data['Cell_PTCL7_484KI67_L']
del data['Cell_PTCL7_484CD11C_L']
del data['Cell_PTCL7_484CD7_L']

'Cell_PTCL7_484KI67_S'
'Cell_PTCL7_484CD11C_S'
'Cell_PTCL7_484CD7_S'

data.rename({'Cell_PTCL7_484KI67_S': 'Cell_PTCL7_484KI67', 
             'Cell_PTCL7_484CD11C_S': 'Cell_PTCL7_484CD11C', 
             'Cell_PTCL7_484CD7_S': 'Cell_PTCL7_484CD7'}, axis=1, inplace=True)

data.to_csv('/Users/aj/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/4-celltype_calling/PTCL7/PTCL7.csv')


#===============================================================================
# Prepare and save PTCL1 dataset
image_path = ['/Users/aj/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/raw_data/whole_sections/PTCL1_450.csv']
adata = histocat_to_object (image_path, remove_string_from_name = "Cell_PTCL1_450", islog=True, drop_markers= ["PERK", "NOS2","BG1","BG2","BG3","ACTIN"])
adata = rescale (adata, gate=None, return_gates=False, failed_markers=['CD15', 'CD56'])
phenotype = pd.read_csv('/Users/aj/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/4-celltype_calling/PTCL1/phenotype_workflow.csv')
adata = phenotype_groups (adata, phenotype, gate = 0.5, label="phenotype")
adata.write('/Users/aj/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/4-celltype_calling/PTCL1/PTCL1.h5ad')
adata.obs['phenotype'].to_csv('/Users/aj/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/4-celltype_calling/PTCL1/pheno.csv')

#===============================================================================
# Prepare and save PTCL2 dataset
data = pd.read_csv('//research.files.med.harvard.edu/ImStor/sorger/data/RareCyte/Connor/Z155_PTCL/whole_section_csv_files/raw/PTCL1_450.csv')
adata = df_to_annobject (data, file_name = "Cell_PTCL1_450", islog=True, drop_markers= ["CellId", "ImageId","PERK", "NOS2","BG1","BG2","BG3","ACTIN"])
adata.write('/Users/aj/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/4-celltype_calling/PTCL2/PTCL2.h5ad')
adata.obs['phenotype'].to_csv('/Users/aj/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/4-celltype_calling/PTCL2/pheno.csv')
#===============================================================================

#===============================================================================
# Prepare and save PTCL3 dataset
data = pd.read_csv('//research.files.med.harvard.edu/ImStor/sorger/data/RareCyte/Connor/Z155_PTCL/whole_section_csv_files/raw/PTCL3_690.csv')
adata = df_to_annobject (data, file_name = "Cell_PTCL3_690", islog=True, drop_markers= ["CellId", "ImageId","PERK", "NOS2","BG1","BG2","BG3","ACTIN"])
adata = rescale (adata, gate=None, return_gates=False)
adata.write('C:/Users/ajn16/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/4-celltype_calling/PTCL3/PTCL3.h5ad')
adata.obs['phenotype'].to_csv('/Users/aj/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/4-celltype_calling/PTCL3/pheno.csv')
#===============================================================================

#===============================================================================
# Prepare and save PTCL4 dataset
data = pd.read_csv('//research.files.med.harvard.edu/ImStor/sorger/data/RareCyte/Connor/Z155_PTCL/whole_section_csv_files/raw/PTCL4_713.csv')
adata = df_to_annobject (data, file_name = "Cell_PTCL4_713", islog=True, drop_markers= ["CellId", "ImageId","PERK", "NOS2","BG1","BG2","BG3","ACTIN"])
adata = rescale (adata, gate=None, return_gates=False)
adata.write('C:/Users/ajn16/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/4-celltype_calling/PTCL4/PTCL4.h5ad')
adata.obs['phenotype'].to_csv('/Users/aj/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/4-celltype_calling/PTCL4/pheno.csv')
#===============================================================================

#===============================================================================
# Prepare and save PTCL5 dataset
data = pd.read_csv('//research.files.med.harvard.edu/ImStor/sorger/data/RareCyte/Connor/Z155_PTCL/whole_section_csv_files/raw/PTCL5_506.csv')
adata = df_to_annobject (data, file_name = "Cell_PTCL5_506", islog=True, drop_markers= ["CellId", "ImageId","PERK", "NOS2","BG1","BG2","BG3","ACTIN"])
adata = rescale (adata, gate=None, return_gates=False)
adata.write('C:/Users/ajn16/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/4-celltype_calling/PTCL5/PTCL5.h5ad')
adata.obs['phenotype'].to_csv('/Users/aj/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/4-celltype_calling/PTCL5/pheno.csv')
#===============================================================================

#===============================================================================
# Prepare and save PTCL6 dataset
image_path = ['/Users/aj/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/raw_data/whole_sections/PTCL6_306.csv']
adata = histocat_to_object (image_path, remove_string_from_name = "Cell_PTCL6_306", islog=True, drop_markers= ["PERK", "NOS2","BG1","BG2","BG3","ACTIN"])
adata = rescale (adata, gate=pd.DataFrame(['CD45', 6.8]).T, return_gates=False, failed_markers=['CD15','CD21','CD25'])
phenotype = pd.read_csv('/Users/aj/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/4-celltype_calling/PTCL6/phenotype_workflow.csv')
adata = phenotype_groups (adata, phenotype, gate = 0.5, label="phenotype")
adata.write('/Users/aj/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/4-celltype_calling/PTCL6/PTCL6.h5ad')
adata.obs['phenotype'].to_csv('/Users/aj/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/4-celltype_calling/PTCL1/pheno.csv')
#===============================================================================

#===============================================================================
# Prepare and save PTCL7 dataset
image_path = ['/Users/aj/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/raw_data/whole_sections/PTCL7.csv']
adata = histocat_to_object (image_path, remove_string_from_name = "Cell_PTCL7_484", islog=True, drop_markers= ["PERK", "NOS2","BG1","BG2","BG3","ACTIN"])
adata = rescale (adata, gate=None, return_gates=False, failed_markers=['CD15','CD21','CD25'])
phenotype = pd.read_csv('/Users/aj/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/4-celltype_calling/PTCL7/phenotype_workflow.csv')
adata = phenotype_groups (adata, phenotype, gate = 0.5, label="phenotype")
adata.write('/Users/aj/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/4-celltype_calling/PTCL7/PTCL7.h5ad')
adata.obs['phenotype'].to_csv('/Users/aj/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/4-celltype_calling/PTCL1/pheno.csv')
#===============================================================================

#===============================================================================
# Prepare and save PTCL8 dataset
data = pd.read_csv('//research.files.med.harvard.edu/ImStor/sorger/data/RareCyte/Connor/Z155_PTCL/whole_section_csv_files/raw/PTCL8_360.csv')
adata = df_to_annobject (data, file_name = "Cell_PTCL8_360", islog=True, drop_markers= ["CellId", "ImageId","PERK", "NOS2","BG1","BG2","BG3","ACTIN"])
adata = rescale (adata, gate=None, return_gates=False)
adata.write('C:/Users/ajn16/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/4-celltype_calling/PTCL8/PTCL8.h5ad')
adata.obs['phenotype'].to_csv('/Users/aj/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/4-celltype_calling/PTCL8/pheno.csv')
#===============================================================================

#===============================================================================
# Prepare and save PTCL9 dataset
data = pd.read_csv('//research.files.med.harvard.edu/ImStor/sorger/data/RareCyte/Connor/Z155_PTCL/whole_section_csv_files/raw/PTCL9_320.csv')
adata = df_to_annobject (data, file_name = "Cell_PTCL9_320", islog=True, drop_markers= ["CellId", "ImageId","PERK", "NOS2","BG1","BG2","BG3","ACTIN"])
adata = rescale (adata, gate=None, return_gates=False)
adata.write('C:/Users/ajn16/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/4-celltype_calling/PTCL9/PTCL9.h5ad')
adata.obs['phenotype'].to_csv('/Users/aj/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/4-celltype_calling/PTCL9/pheno.csv')
#===============================================================================

#===============================================================================
# Prepare and save PTCL10 dataset
data = pd.read_csv('//research.files.med.harvard.edu/ImStor/sorger/data/RareCyte/Connor/Z155_PTCL/whole_section_csv_files/raw/Ton_192.csv')
adata = df_to_annobject (data, file_name = "Cell_Ton_192", islog=True, drop_markers= ["CellId", "ImageId","PERK", "NOS2","BG1","BG2","BG3","ACTIN"])
adata = rescale (adata, gate=None, return_gates=False)
adata.write('C:/Users/ajn16/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/4-celltype_calling/PTCL10/PTCL10.h5ad')
adata.obs['phenotype'].to_csv('/Users/aj/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/4-celltype_calling/PTCL10/pheno.csv')
#===============================================================================

