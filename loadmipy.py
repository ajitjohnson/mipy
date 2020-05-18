#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 08:55:34 2020
@author: Ajit Johnson Nirmal
Load MiPy
"""
import sys

#  Import SiPy package components
sys.path.insert(1, '/Users/aj/Dropbox (Partners HealthCare)/packages/mipy')
#sys.path.insert(1, 'C:/Users/ajn16/Dropbox (Partners HealthCare)/packages/mipy')
from pp.add_imageid import add_imageid
from pp.mcmicro_to_object import mcmicro_to_object
from pp.rescale import rescale
from tl.phenotype_groups import phenotype_groups
from tl.phenotype_subcluster import phenotype_subcluster
from tl.rename_clusters import rename_clusters
from tl.add_meta import add_meta
from tl.community_clustering import community_clustering
from tl.cell_proportion import cell_proportion
from tl.cell_abundance_corr import *
from tl.knn import knn
from tl.community_detection import community_detection
from tl.community_clustering import community_clustering
from tl.community_umap import community_umap
from tl.nn_interaction import nn_interaction
from tl.proximal_density_score import proximal_density_score
from tl.adata_concatenate import adata_concatenate


from pl.marker_expression_plot import marker_expression_plot
from pl.plotly import plotly
from pl.percent_plot import percent_plot
from pl.nn_interaction_heatmap import nn_interaction_heatmap
from pl.napari import *
from pl.pca_het import pca_het
from pl.gmm_dist_plot import gmm_dist_plot
from pl.gate_finder import gate_finder
from export.save_files import *

