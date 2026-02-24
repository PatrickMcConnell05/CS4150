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

smallest = 100000000000
largest = 0

for col in df.columns[3:]:
    values = df[col] #sets values to be a single column
    count = (values > 0).sum() #gets the number of 1s in a col
    if(count < smallest):
        smallest = count
    if(count > largest):
        largest = count

# print("Smallest Number of Windows:", smallest)
# print("Largest Number of Windows:", largest)

group_size = (largest - smallest) / 5
# print(group_size)


#6.) work from act2.py --meant to pair names and values up

np_names = df.columns.values[3:]

df2 = df.loc[:, df.columns.intersection(np_names)].copy() #this creates the dataframe (df3) with the nps that had at least 1 window in the hist1 section
# print(df3)

count = 0

arr_windows_per_np = []


for col in df2.columns: #goes through all of the columns
    values = df2[col] #assigns a col to values
    count += (values > 0).sum() #sums the values
    arr_windows_per_np.append(count) #adds the 
    count = 0

    
pairing = dict(zip(np_names, arr_windows_per_np)) #key value pairing of the np names and the num of windows present
sorted_pairings = dict(sorted(pairing.items(), key=lambda item: item[1]))

strongly_apical = []
somewhat_apical = []
neither = []
somewhat_equatorial = []
strong_equatorial = []

pairs = list(sorted_pairings.items())

for np, frequency in pairs:
    if frequency <= group_size:
        strongly_apical.append(((np, frequency)))
    elif frequency <= (group_size * 2):
        somewhat_apical.append((np,frequency))
    elif frequency <= (group_size * 3):
        neither.append((np,frequency))
    elif frequency <= (group_size * 4):
        somewhat_equatorial.append((np,frequency))
    # elif frequency <= (group_size * 5):
    #     strong_equatorial.append((np,frequency))
    else:
        strong_equatorial.append((np,frequency))
        

for index, cluster in enumerate(best_clusters):
    
    strongly_apical_count = 0
    somewhat_apical_count = 0
    neither_count = 0
    somewhat_equatorial_count = 0
    strong_equatorial_count = 0

    for np_name in cluster:
        
        for i in range(len(cluster)):
            # if cluster[np_name] in strongly_apical:
            #     strongly_apical_count += 1
            # elif cluster[np_name] in somewhat_apical:
            #     somewhat_apical += 1
            # elif cluster[np_name] in neither:
            #     neither_count += 1
            # elif cluster[np_name] in somewhat_equatorial:
            #     somewhat_equatorial_count += 1
            # elif cluster[np_name] in strong_equatorial:
            #     strong_equatorial_count += 1
            if np_name in dict(strongly_apical).keys():
                strongly_apical_count += 1
            elif np_name in dict(somewhat_apical).keys():
                somewhat_apical_count += 1
            elif np_name in dict(neither).keys():
                neither_count += 1
            elif np_name in dict(somewhat_equatorial).keys():
                somewhat_equatorial_count += 1
            elif np_name in dict(strong_equatorial).keys():
                strong_equatorial_count += 1
    
    strongly_aplical_percentage = strongly_apical_count / len(cluster)
    somewhat_apical_percentage = somewhat_apical_count / len(cluster)
    neither_percentage = neither_count / len(cluster)
    somewhat_equatorial_percentage = somewhat_equatorial_count / len(cluster)
    strong_equatorial_percentage = strong_equatorial_count / len(cluster)
    
    plt.figure(figsize=(10,6))
    plt.bar("Strongly Apical", strongly_aplical_percentage)
    plt.bar("Somewhat Apical", somewhat_apical_percentage)
    plt.bar("Neither", neither_percentage)
    plt.bar("Somewhat Equatorial", somewhat_equatorial_percentage)
    plt.bar("Strongly Equatorial", strong_equatorial_percentage)
    plt.title(f"Radial Position Distribution for Cluster {index}")
    plt.ylabel("Percentage of NPs")
    plt.savefig(f"./bar_graphs/cluster_{index}_radial_distribution.png")
    plt.close()
        












############## --- ACTIVITY 3 --- ###############





