import sys
sys.path.append("..")

from functions import *
import numpy as np
import pandas as pd


data_path = "../data/GSE64881_segmentation_at_30000bp.passqc.multibam.txt"
feat_path = "../data/Hist1_region_features.csv"
df, np_cols = load_and_filter(data_path)
feat = load_features(feat_path) #loads the features for the nodes, which will be used to calculate the percentage of nodes in the community that contain a hist1 gene and the percentage of nodes in the community that contain a LAD
norm_linkage_matrix = compute_normalized_linkage_matrix(df, np_cols) #computes norm link matrix

#explain vals very briefly:
#vals in an array of the upper triangle values of the normalized linkage matrix, sorted in ascending order
vals = norm_linkage_matrix.values[np.triu_indices_from(norm_linkage_matrix.values,k=1)]
vals.sort()
Q3 = np.percentile(vals,75) #gets the value to seperate the 75th percentile

adj_matrix = create_adj_matrix(norm_linkage_matrix, Q3)

dc = degree_centrality(adj_matrix)
top_five_dc = top_five_degree_centrality(dc) #gets the top 5 nodes with highest degree centrality, with their names and values
# print("Top 5 nodes with highest degree centrality:")
# print(top_five_dc)




clusters = cluster_by_top_degree_centrality(adj_matrix, norm_linkage_matrix, top_five_dc.index)
print_community_information(clusters, feat)

community_visualization(clusters, adj_matrix, dc)
community_heatmap(clusters, norm_linkage_matrix)

