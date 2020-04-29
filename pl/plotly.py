#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 17:10:46 2020
@author: Ajit Johnson Nirmal
Plotly for visualising Phenotypes
"""

# Vizualising using plotly
import pandas as pd
import plotly.express as px
import plotly.io as pio
pio.renderers.default = 'browser'


def plotly (adata,phenotype,image_id=None,x='X_position',y='Y_position',size=2):
    if image_id is not None:
        adata = adata[adata.obs['ImageId'] == image_id]    
    data = pd.DataFrame({'x':adata.obs[x], 'y':adata.obs[y],'col': adata.obs[phenotype]})
    fig = px.scatter(data, x="x", y="y", color="col")
    fig.update_traces(marker=dict(size=size),selector=dict(mode='markers'))
    fig.update_yaxes(autorange="reversed")
    return fig