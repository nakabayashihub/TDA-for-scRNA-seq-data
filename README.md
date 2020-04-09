# TDA-for-scRNA-seq-data
Single cell RNA-seq data analysis using Topological Data Analysis (TDA)  
Topology is one of a branch of mathematics that deals with the mathematical property of shapes'.  Topological Data Analysis (TDA) is recently focused as a promising method provides a insight into the structure underling a data cloud in high dimensional space by the technics from topology.  TDA has found many applications in many fields including biology, visual recognition, natural language processing, and so on.  Here scRNA-seq data is analyzed using TDA.
## Dataset
[GSE67310](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67310)  
Dataset of scRNA-seq for 405 induced neuronal (iN) cells from Mouse Embryonic Fibroblasts (MEFs) at five time points.  
MEFs are reprogramed by over expression of Ascl1 gene.  Continuous reprograming path is reconstructed through transriptome similarity.
|time point|#cells|
|:------:|:----:|
|day 0|73|
|day 2|128|
|day 5|55|
|day 20|33|
|day 22|116|

## python code  
First, Vietris-Rips filtration is performed to compute persistent homology.  A network composed of cells is constructed from 1-simplex that is to say linkage between data points.  Obtained network is represented using netwrokx library as undirected graph.  Circles represent cells.  Time points or types of cells are indeicated by color.
~~~python
#Load libraries
import numpy as np
import pandas as pd
import gudhi as gd
import networkx as nx
#Reading data
x = pd.read_csv('GSE67310_iN_data_log2FPKM_annotated.txt', delimiter = '\t')
#Triming data
y = x.drop('cell_name', axis = 1)
y = y.drop('assignment', axis = 1)
y = y.drop('log_tauGFP_intensity', axis = 1)
y = y.drop('experiment', axis = 1)
y = y.drop('time_point', axis = 1)
y.index = x.cell_name
#Creating color table by day
day_color = pd.DataFrame()
for i in range(len(x)):
  if x.time_point[i] == 0:
    day_color[y.index[i]] = 'red'
  elif x.time_point[i] == 2:
    day_color[y.index[i]] = 'yellow'
  elif x.time_point[i] == 5:
    day_color[y.index[i]] = 'green'
  elif x.time_point[i] == 20:
    day_color[y.index[i]] = 'purple'
  else:
     day_color[y.index[i]] = 'blue'
#Creating color table by cell type
type_color = pd.Series()
for i in range(len(x)):
  if x.assignment[i] == 'MEF':
    type_color[y.index[i]] = 'red'
  elif x.assignment[i] == 'd2_induced':
    type_color[y.index[i]] = 'yellow'
  elif x.assignment[i] = 'd2_intermediate':
    type_color[y.index[i]] = 'orange'
  elif x.assignment[i] == 'd5_earlyiN':
    type_color[y.index[i]] = 'skyblue'
  elif x.assignment[i] == 'd5_earlyMyocyte':
    type_color[y.index[i]] = 'lightgeen'
  elif x.assignment[i] == 'd5_intermediate':
    type_color[y.index[i]] = 'brown'
  elif x.assignment[i] == 'd5_failedReprog':
    type_color[y.index[i]] = 'gray'
  elif x.assignment[i] == 'd22_failedReprog':
    type_color[y.index[i]] = 'black'
  elif x.assignment[i] == 'Neuron':
    type_color[y.index[i]] = 'blue'
  elif x.assignment[i] == 'Myocyte':
    type_color[y.index[i]] = 'green'
  else:
    type_color[y.index[i]] = 'white'
#Computing Vietris-Rips complex
rips = gd.RipsComplex(y.values, max_edge_length = 250)
#Computing simplex tree
simplex_tree = rips.create_simplex_tree(max_dimension = 2)
#Computing skeleton
skeleton = simplex_tree.get_skeleton(2)
#Getting persistence diagram
diag = simplex_tree.persistence()
#Plotting persistence diagram
gd.plot_persistence_diagram(diag)
#Plotting persistence density
gd.plot_persistence_density(diag)
#Constructing network
g = nx.Graph()
for i in range(len(skeleton)):
  if len(skeleton[i][0]) == 2:
    g.add_edge(y.index[skeleton[i][0][0]], y.index[skeleton[i][0][1]])
layout = nx.kamada_kawai_layout(g)
nx.draw_networkx_nodes(g, layout, linewidths = 0.2, edgecolors = 'black', node_size = 20, node_color = day_color[list(g.nodes())].values)
nx.draw_networkx_edges(g, layout, width = 0.2, edge_color = 'gray')
nx.draw_networkx_nodes(g, layout, linewidths = 0.2, edgecolors = 'black', node_size = 20, node_color = type_color[list(g.nodes()].values)
nx.draw_networkx_edges)g, layout, width = 0.2, edge_color = 'gray')
~~~
![persistence_diagram](iNeuron_Persistence_diagram.tiff)
![persistence_density](iNeuron_Persistence_density.tiff)
![network_day](iNeuron_Network_day_samll.tiff)
![network_type](iNeuron_Network_type_small.tiff)
