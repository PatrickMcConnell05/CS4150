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

best_centers, best_clusters, best_var = best_of_n_runs(norm_distance_matrix)

print("Final Centers:", best_centers)
print("Lowest Within Cluster Variance:", best_var)


############## --- ACTIVITY 2 --- ###############
feat_table_path = "./data/Hist1_region_features.csv"
feat_table = pd.read_csv(feat_table_path)

for cluster in best_clusters:
    for np in cluster:
        for window in df.loc[np]:
            if feat_table.loc[window][df] == 1 and np[window] == 1:
                count += 1




############## --- ACTIVITY 2 --- ###############





