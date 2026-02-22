import pandas as pd
import numpy as np
import random
import os
from functions import *

data_path = "./data/GSE64881_segmentation_at_30000bp.passqc.multibam.txt"

#data loading and preprocessing
df, np_cols = load_and_filter(data_path)
jaccard_matrix, dist_matrix = compute_jaccard(df, np_cols)
norm_jaccard_matrix, norm_distance_matrix = compute_normalized_jaccard(df, np_cols)


three_points = list(np.random.choice(np_cols, size=3, replace=False))
print("Original Three Points:",three_points)

one_run = run_kmedoids(norm_distance_matrix, three_points)

print("Variation for original points:", one_run["init_var"])
print("Final Total Within Cluster Variation:", one_run["var"])
print("Final Centers:", one_run["centers"])
print("Cluster 1:", len(one_run["clusters"][0]))
print("Cluster 2:", len(one_run["clusters"][1]))
print("Cluster 3:", len(one_run["clusters"][2]))


best_centers, best_clusters, best_var = best_of_n_runs(norm_distance_matrix)

print("Final Centers:", best_centers)
print("Lowest Within Cluster Variance:", best_var)