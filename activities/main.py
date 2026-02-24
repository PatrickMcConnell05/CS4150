import pandas as pd
import numpy as np
import random
import os
import matplotlib.pyplot as plt
from functions import *

data_path = "./data/GSE64881_segmentation_at_30000bp.passqc.multibam.txt"
hist1_feat_path = "./data/Hist1_region_features.csv"

#data loading and preprocessing
df, np_cols = load_and_filter(data_path)
jaccard_matrix, dist_matrix = compute_jaccard(df, np_cols)
norm_jaccard_matrix, norm_distance_matrix = compute_normalized_jaccard(df, np_cols)


############## --- ACTIVITY 1 --- ###############
three_points = list(np.random.choice(np_cols, size=3, replace=False))
print("Original Three Points:",three_points)

one_run = run_kmedoids(norm_distance_matrix, three_points)

print("Variation for original points:", one_run["init_var"])

best_centers, best_clusters, best_var = best_of_n_runs(norm_distance_matrix)

print("Final Centers:", best_centers)
print("Lowest Within Cluster Variance:", best_var)

############## --- ACTIVITY 1 --- ###############



############## --- ACTIVITY 2 --- ###############
feat_table_path = "./data/Hist1_region_features.csv"
feat = load_features(feat_table_path) #loads in the feature table, plus a chrom, start, and stop columns at the end


#df with: cluster, NP, Hist1_pct, and LAD_pct columns
results_df = compute_feature_percentages(df, feat, best_clusters) 
print(results_df)



plot_feature_boxplots(
    results_df,
    "Hist1_pct",
    "% of NP windows with Hist1 genes",
    "Hist1 Feature by Cluster",
    "./heatmaps/hist1_boxplot.png"
)

plot_feature_boxplots(
    results_df,
    "LAD_pct",
    "% of NP windows with LADs",
    "LAD Feature by Cluster",
    "./heatmaps/lad_boxplot.png"
)



############## --- ACTIVITY 2 --- ###############






