import sys
sys.path.append("..")

from functions import *
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', None)


#data loading and preprocessing
data_path = "../data/GSE64881_segmentation_at_30000bp.passqc.multibam.txt"
df, np_cols = load_and_filter(data_path)
norm_linkage_matrix = compute_normalized_linkage_matrix(df, np_cols) #computes norm link matrix




vals = norm_linkage_matrix.values[np.triu_indices_from(norm_linkage_matrix.values,k=1)]
vals.sort()

Q3 = np.percentile(vals,75) #gets the value to seperate the 75th percentile

# print(Q3)
# print(np_cols)

adj_matrix = create_adj_matrix(norm_linkage_matrix, Q3)
# print(adj_matrix)

dc = degree_centrality(adj_matrix)
# print_degree_centrality_stats(dc)
network_graph(adj_matrix, "output.png")