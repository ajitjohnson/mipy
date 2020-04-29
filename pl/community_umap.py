#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 23:24:38 2020
@author: Ajit Johnson Nirmal
UMAP of the cell communities
"""

import anndata as ad
import scanpy as sc
import pandas as pd


def community_umap (adata,community_matrix='nn_communities',color=None,subset_image=None,
                    unique_id='ImageId',n_neighbors=30):
    
    # Subset a particular image if needed
    if subset_image is not None:
        bdata = adata[adata.obs[unique_id] == subset_image]
    else:
        bdata = adata
    
    d = bdata.uns[community_matrix]
    d = d.loc[bdata.obs.index].fillna(0)  
    d = d[(d.T != 0).any()] # Drop all rows with zero   
    samples = d.index  # Get the cells that survived
    cdata = bdata[samples] # Do this so that the .obs of this can be used below
    
    # Create a new ann data object
    bdata = ad.AnnData(d)
    bdata.obs = pd.DataFrame(cdata.obs)
    
    # Convert umap coloring to category
    if color is not None:
        if color in list(bdata.obs.columns):
            bdata.obs[color] = bdata.obs[color].astype('category')
    
    # Build the UMAP
    sc.tl.pca(bdata)
    sc.pp.neighbors(bdata, n_neighbors=n_neighbors)
    sc.tl.umap(bdata)
    
    # Add the UMAP embedding to adata
    adata.obsm['X_umap'] = bdata.obsm['X_umap']
    
    # Return
    return adata




adata.obsm['t'] = d['B cells']
    
    
    sc.pl.umap(bdata, color=color, size = size, use_raw=False, color_map=color_map,legend_loc=legend_loc) #legend_loc='on data'
    


    sc.tl.pca(adata)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors)
    sc.tl.umap(adata)
    sc.pl.umap(adata, color=color, size = size, use_raw=False, color_map=color_map,legend_loc=legend_loc) #legend_loc='on data'
    


    sc.pl.umap(adata, color=[], size = size, use_raw=False, color_map=color_map,legend_loc=legend_loc) #legend_loc='on data'


adata.obs['tumor'] = d['Tumor CD30+']


     bdata.var
     
     adata.var
        


bdata.va = bdata.X
    

            



#sc.tl.diffmap(bdata, n_comps=5)
    
        
sc.settings.n_jobs=4


sns.set_style("white")

sc.pl.umap(bdata, color=['original_tma_id'], size = 20, use_raw=False, color_map= 'vlag') #legend_loc='on data'
sc.pl.umap(bdata, color=['ImageId'], size = 10, use_raw=False, color_map= 'vlag') #legend_loc='on data'
sc.pl.umap(bdata, color=['Tumor CD30+'], size = 10, use_raw=False, color_map= 'vlag') #legend_loc='on data'
sc.pl.umap(bdata, color=['M2 Macrophages'], size = 10, use_raw=False, color_map= 'vlag') #legend_loc='on data'


sc.pl.diffmap(bdata, color=['ImageId'], size = 30, use_raw=False, color_map= 'vlag') #legend_loc='on data'

# Cluster the data
sc.tl.leiden(bdata)
sc.pl.umap(bdata, color=['ImageId','leiden'], size = 10, use_raw=False, color_map= 'vlag',legend_loc='on data') #legend_loc='on data'

# Investigate cluster
cluster = 3
bdata [bdata.obs['leiden'] == '3'].obs['ImageId'].value_counts()


adata.obs['leiden'] = bdata.obs['leiden']
cdata = adata [adata.obs['ImageId'] == 215]
m = cdata.var.index
sc.tl.dendrogram(cdata,groupby='leiden')
sc.pl.matrixplot(bdata, var_names= m, groupby='leiden', dendrogram=True, use_raw=False, cmap='viridis_r',vmin=0,vmax=1, swap_axes=True)


cdata = bdata [bdata.obs['ImageId'] == 215]
m = cdata.var.index
sc.tl.dendrogram(cdata,groupby='leiden')
sc.pl.matrixplot(cdata, var_names= m, groupby='leiden', dendrogram=True, use_raw=False, cmap='viridis_r',vmin=0,vmax=0.8, swap_axes=True)
sc.pl.heatmap(cdata, var_names=m, groupby='leiden', swap_axes=True, use_raw=False, log=False, cmap= 'vlag')

sc.tl.rank_genes_groups(bdata, 'leiden', method='t-test')
sc.pl.rank_genes_groups(bdata, n_genes=5, sharey=False)

x = pd.DataFrame(bdata.uns['rank_genes_groups']['names']).head(5)


sc.tl.pca(bdata)
sc.pl.pca(bdata, color='Tumor',color_map= 'vlag')

bdata.obs[phenotype].unique()