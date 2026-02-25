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

BEST_CLUST_CENTERS = ['F10C3', 'F16F4', 'F7F3']


############## --- ACTIVITY 1 --- ###############
# three_points = list(np.random.choice(np_cols, size=3, replace=False))
# print("Original Three Points:",three_points)

# one_run = run_kmedoids(norm_distance_matrix, three_points)

# print("Variation for original points:", one_run["init_var"])

# best_centers, best_clusters, best_var = best_of_n_runs(norm_distance_matrix)

# print("Final Centers:", best_centers)
# print("Lowest Within Cluster Variance:", best_var)


### shortended version of act 1
best_clusters = assign_to_k_clusters(norm_distance_matrix, BEST_CLUST_CENTERS)
# print("best_clusters:", best_clusters)
###

############## --- ACTIVITY 1 --- ###############



############## --- ACTIVITY 2 --- ###############
feat_table_path = "./data/Hist1_region_features.csv"
feat = load_features(feat_table_path) #loads in the feature table, plus a chrom, start, and stop columns at the end


#df with: cluster, NP, Hist1_pct, and LAD_pct columns
results_df = compute_feature_percentages(df, feat, best_clusters) 
# print(results_df)



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



############## --- ACTIVITY 3 --- ###############

# 1) compute windows-per-NP
df2 = df[np_cols].copy()

windows_per_np = {}
for col in df2.columns:
    windows_per_np[col] = int((df2[col] > 0).sum())

smallest = min(windows_per_np.values())
largest  = max(windows_per_np.values())

radial_pos = {}
if largest == smallest:
    radial_pos = {np: 3 for np in windows_per_np}  # or 1, but 3 is more “neutral”
else:
    group_size = (largest - smallest) / 5.0
    for np, freq in windows_per_np.items():
        if freq <= smallest + group_size:
            radial_pos[np] = 1
        elif freq <= smallest + 2*group_size:
            radial_pos[np] = 2
        elif freq <= smallest + 3*group_size:
            radial_pos[np] = 3
        elif freq <= smallest + 4*group_size:
            radial_pos[np] = 4
        else:
            radial_pos[np] = 5

labels = [
    "1 Strongly apical",
    "2 Somewhat apical",
    "3 Neither",
    "4 Somewhat equatorial",
    "5 Strongly equatorial"
]


os.makedirs("./bar_graphs", exist_ok=True) #nsures that the directory exists for saving the bar graphs

# 3) for each cluster, count + percentage + plot
for ci, cluster in enumerate(best_clusters, start=1):
    counts = [0, 0, 0, 0, 0]

    for np_name in cluster:
        rp = radial_pos.get(np_name, None)
        if rp is not None:
            counts[rp - 1] += 1

    total = sum(counts) if sum(counts) > 0 else 1
    perc = [c / total for c in counts]

    # Bar height should be NUMBER of NPs (per assignment)
    plt.figure(figsize=(10, 6))
    plt.bar(labels, counts)
    plt.title(f"Radial Position Counts for Cluster {ci}")
    plt.ylabel("Number of NPs")
    plt.xticks(rotation=20, ha="right")
    plt.tight_layout()
    plt.savefig(f"./bar_graphs/cluster_{ci}_radial_counts.png")
    plt.close()

    # optional: print percentages for your writeup
    print(f"\nCluster {ci} (n={total}) radial %:")
    for lab, p in zip(labels, perc):
        print(f"  {lab}: {p*100:.1f}%")
        












############## --- ACTIVITY 3 --- ###############





