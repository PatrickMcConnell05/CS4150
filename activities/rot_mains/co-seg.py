from functions import *
data_path = "./data/GSE64881_segmentation_at_30000bp.passqc.multibam.txt"

#data loading and preprocessing
df, np_cols = load_and_filter(data_path)

norm_linkage_matrix = compute_normalized_linkage_matrix(df, np_cols) #computes norm link matrix

save_heatmap(
    norm_linkage_matrix, 
    "Normalized Linkage Matrix - Hist1 Region", 
    "./heatmaps/co-segregation/normalized_linkage.png", 
    cbar_label = "Normalized Linkage",
    ticks=[-1, -0.5, 0, 0.5, 1],
    vmin=-1,
    vmax=1,
    xticklabels=False,
    yticklabels=False
)

print(norm_linkage_matrix)
