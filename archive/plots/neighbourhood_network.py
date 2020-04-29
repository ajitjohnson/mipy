#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 17:55:31 2020

@author: aj
"""

import networkx as nx

corr = neighbours_max.corr()

corr = neighbours_max

links = corr.stack().reset_index()
links.columns = ['var1', 'var2','value']
links


links_filtered=links.loc[ (links['value'] > 0.5) & (links['var1'] != links['var2']) ]
links_filtered


G=nx.from_pandas_edgelist(links_filtered, 'var1', 'var2')
G = nx.from_pandas_edgelist(links_filtered,source='var1',target='var2', edge_attr=None, create_using=nx.DiGraph())


nx.draw(G, with_labels=True, node_color='orange', node_size=400, edge_color='black', linewidths=2, font_size=7,pos=nx.fruchterman_reingold_layout(G))




# New neighbourhood analysis
