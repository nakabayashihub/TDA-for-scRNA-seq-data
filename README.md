# TDA-for-scRNA-seq-data
Single cell RNA-seq data analysis using Topological Data Analysis (TDA)
Questions?Comments?Contact
#Dataset
GSE67310
~~~python
#load libraries
import numpy as np
import pandas as pd
import gudhi as gd
import networkx as nx
#reading data
x = pd.read_csv('GSE67310_iN_data_log2FPKM_annotated.txt', delimiter = '\t')
#triming data
y = x.drop('cell_name', axis = 1)
y = x.drop('assignment', axis = 1)
y = x.drop('log_tauGFP_intensity', axis = 1)
y = x.drop('experiment', axis = 1)
y = x.drop('time_point', axis = 1)
#Computing Vietris-Rips complex
rips = gd.RipsComplex(y.values, max_edge_lemgth = 200)
#Computing simplex tree
simplex_tree = rips.create_simplex_tree(max_dimension = 2)
#Computing skeleton
skeleton = simplex_tree.get_skeleton(2)
#Constructing netowrk
g = nx.Graph()
for i in range(len(skeleton)):
  if len(skeleton[i][0]) == 2:
    g.add_edges(skeleton[i][0][0], skeleton[i][0][1])
layout = nx.kamada_kawai_layout(g)
nx.draw_networkx_node(g,layout,lineidths=0.2, edgecolors='black', node_size=20, node_color
~~~
