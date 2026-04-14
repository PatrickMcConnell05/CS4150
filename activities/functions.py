import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os







#takes in center points and the distances to centers and assigns clusters 
def assign_to_k_clusters(dist_to_centers: pd.DataFrame, centers: list[str]):
    
    k = len(centers) #will always just be three if we select 3 center points
    clusters = [[] for _ in range(k)] #creates a list of lists of size k (will just be the three clusters)
    
    for np_name in dist_to_centers.index:
        distances = []

        for i in range(len(centers)):
            d = dist_to_centers.loc[np_name, centers[i]] #distance of to each np center in each row
            distances.append(d)
        
        min_dist = min(distances)
        min_index = distances.index(min_dist) #gets the index of the np with the smallest distance to a given np
        
        clusters[min_index].append(np_name) #adds the np at the index of the cluster with the smallest distance
        
    return clusters






#finds and returns the np in a cluster that has the minimum disssimilarity
def find_center(cluster: list[str], norm_dist_matrix: pd.DataFrame) -> str:
    if len(cluster) == 0:
        return None

    best_np = cluster[0]
    best_avg_diss = float("inf")
    
    for np1 in cluster: #loops through all of the nps
        total = 0.0
        for np2 in cluster: #again loops through all nps
            total += norm_dist_matrix.loc[np1,np2] #gets distance and adds to ttoal
        avg_diss = total / len(cluster)

        if avg_diss < best_avg_diss: #if a new lowest is found update
            best_avg_diss = avg_diss
            best_np = np1
            
    return best_np










#sums up the distances between all points in a cluster to the center, and sums for all clusters
def within_cluster_var(clusters: list[list[str]], centers: list[str], norm_dist_matrix: pd.DataFrame):
    total = 0.0 #the total distances from centers, within cluster variation
    
    for i, cluster in enumerate(clusters): #goes through all of the clusters
        if len(cluster) == 0:
            continue

        cent = centers[i]
        total += norm_dist_matrix.loc[cluster, cent].sum() #all of the nps at the cluster center column summed

    return total










#Reclusters until cluters are unchanging or specified number of iterations
def cluster_medoids(norm_dist_matrix: pd.DataFrame, original_centers: list[str], max_iter: int=1000):
    centers = list(original_centers)
    
    for it in range(max_iter):
        dist_to_centers = norm_dist_matrix[centers] #163x3
        clusters = assign_to_k_clusters(dist_to_centers, centers)
        
        
        #new centers
        new_centers = []
        for i, cluster in enumerate(clusters): #computes the new medoid for each cluster
            c = find_center(cluster, norm_dist_matrix) #c should have min. avg. dissimilarity
            if c is None: #handles empty cluster, keeps previous center 
                c = centers[i]
            new_centers.append(c)
        
        new_clusters = assign_to_k_clusters(norm_dist_matrix[new_centers], new_centers) #reclusters based on new centers

        # wc_var = within_cluster_var(new_clusters, new_centers, norm_dist_matrix) #gets the within cluster var for the new clusters
        
        #output --removed for the FEATURE SELECTION - ACTIVITY1
        # print("\nIteration:", it+1)
        # print("Variation:",wc_var)
        # print(new_centers, '\n')

        if set(new_centers) == set(centers):
            centers = new_centers
            clusters = new_clusters
            break
        
        
        
        centers = new_centers
        clusters = new_clusters
    
    return centers, clusters #returns the updated centers and clusters   










def sample_valid_centers(norm_dist_matrix: pd.DataFrame, k: int = 3, max_tries: int = 10000):
    points = list(norm_dist_matrix.index)

    for _ in range(max_tries):
        centers = list(np.random.choice(points, size=k, replace=False))
        clusters = assign_to_k_clusters(norm_dist_matrix[centers], centers)

        if all(len(c) > 0 for c in clusters):
            return centers, clusters

    raise RuntimeError("Couldn't find an initialization with no empty clusters.")



#filters the data and loads it in, returning the df and the list of NP columns, removing any NP columns that have no windows present (all zeros)
#and removing the first three columns (chrom, start, stop) and keeping just the hist1 region (chr13 range we used)
def load_and_filter(path):
    df = pd.read_csv(path, sep="\t")

    df = df[(df["chrom"] == "chr13") &
            (df["start"] >= 21690000) &
            (df["stop"] <= 24120000)].copy()

    # drops NP columns that have no windows present (all zeros)
    for col in df.columns[3:]:
        if df[col].sum() == 0:
            df = df.drop(columns=col)

    np_cols = list(df.columns[3:])
    return df, np_cols










def compute_jaccard(df, np_cols):
    # regular jaccard similarity matrix (m11 / (m01+m10+m11))
    n = len(np_cols)

    jaccard_matrix = pd.DataFrame(
        np.eye(n),
        index=np_cols,
        columns=np_cols
    )

    for a in range(n):
        for b in range(a, n):
            m_01 = m_10 = m_11 = 0

            colA = df[np_cols[a]]
            colB = df[np_cols[b]]

            for x in range(len(df)):
                if colA.iat[x] == 1 and colB.iat[x] == 1:
                    m_11 += 1
                elif colA.iat[x] == 1 and colB.iat[x] == 0:
                    m_10 += 1
                elif colA.iat[x] == 0 and colB.iat[x] == 1:
                    m_01 += 1

            denom = m_01 + m_10 + m_11
            jaccard = (m_11 / denom) if denom > 0 else 0

            jaccard_matrix.iat[a, b] = jaccard
            jaccard_matrix.iat[b, a] = jaccard

    dist_matrix = 1 - jaccard_matrix
    return jaccard_matrix, dist_matrix










def compute_normalized_jaccard(df, np_cols):
    # normalized jaccard similarity matrix (m11 / min(#ones in A, #ones in B))
    n = len(np_cols)

    norm_jaccard_matrix = pd.DataFrame(
        np.eye(n),
        index=np_cols,
        columns=np_cols
    )

    for a in range(n):
        for b in range(a, n):
            num_ones_a = 0
            num_ones_b = 0
            m_11 = 0

            colA = df[np_cols[a]]
            colB = df[np_cols[b]]

            for x in range(len(df)):
                if colA.iat[x] == 1 and colB.iat[x] == 1:
                    m_11 += 1
                if colA.iat[x] == 1:
                    num_ones_a += 1
                if colB.iat[x] == 1:
                    num_ones_b += 1

            denom = min(num_ones_a, num_ones_b)
            jaccard_sim = (m_11 / denom) if denom > 0 else 0

            norm_jaccard_matrix.iat[a, b] = jaccard_sim
            norm_jaccard_matrix.iat[b, a] = jaccard_sim

    norm_dist_matrix = 1 - norm_jaccard_matrix
    return norm_jaccard_matrix, norm_dist_matrix
    
    
    
    
    
    
    
def save_heatmap(matrix, title, out_path, cbar_label, ticks=None,
                 figsize=(12, 10), vmin=0, vmax=1,
                 xticklabels=False, yticklabels=False):
    import matplotlib.pyplot as plt
    import seaborn as sns

    plt.figure(figsize=figsize)
    sns.heatmap(
        matrix,
        cmap="plasma",
        square=True,
        vmin=vmin,
        vmax=vmax,
        xticklabels=xticklabels,
        yticklabels=yticklabels,
        cbar_kws={"label": cbar_label, **({"ticks": ticks} if ticks else {})}
    )
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    
    
    
    
    
    
    
    
    
    
#does one full clustering run starting from given centers
def run_kmedoids(norm_dist_matrix, centers, max_iter=1000):
    
    #initial clusters and variation
    clusters = assign_to_k_clusters(norm_dist_matrix[centers], centers)
    init_var = within_cluster_var(clusters, centers, norm_dist_matrix)

    new_centers, new_clusters = cluster_medoids(norm_dist_matrix, centers, max_iter=max_iter)
    final_var = within_cluster_var(new_clusters, new_centers, norm_dist_matrix)

    return {
        "init_centers": centers,
        "init_clusters": clusters,
        "init_var": init_var,
        "centers": new_centers,
        "clusters": new_clusters,
        "var": final_var
    }
    
    
    
    
    
    
    
    
    
    
    
#tries n number of random initalizations to find clustering croups with lowest within cluster variance
def best_of_n_runs(norm_dist_matrix, n_runs=1000, k=3):
    best_centers = None
    best_clusters = None
    best_var = float("inf")

    for i in range(n_runs):
        print("Iteration:", i + 1)
        centers, clusters = sample_valid_centers(norm_dist_matrix, k=k) #grabbing 3 random points
        centers, clusters = cluster_medoids(norm_dist_matrix, centers) #clusters around the points

        wc_var = within_cluster_var(clusters, centers, norm_dist_matrix) #calculates within cluster variance

        if wc_var < best_var: #checks if a better variance has been found
            best_var = wc_var
            best_centers = centers
            best_clusters = clusters
            print("Variance:",best_var)
            print("Centers:",best_centers)

    return best_centers, best_clusters, best_var



######## ACTIVITY 2 FUNCTIONS ###########

#grabs the feat table cvs and returns
def load_features(feat_path):
    feat = pd.read_csv(feat_path) #takes in the feature table

    #puts the chrom name in the first index
    name_split = feat["name"].str.split(":", expand=True) 

    #splits the 21690000 and 21720000 parts of "21690000-21720000" apart
    range_split = name_split[1].str.split("-", expand=True)

    feat["chrom"] = name_split[0]
    feat["start"] = range_split[0].astype(int)
    feat["stop"]  = range_split[1].astype(int)

    return feat

#computes the Hist1 and LAD percentages for each NP in each cluster, returns a df with cluster, NP, Hist1_pct, and LAD_pct columns
def compute_feature_percentages(df, feat, clusters):
#computes Hist1 and LAD percentages per cluster
    #df = segregation table (with NP columns)
    #feat = feature table aligned to df (must have Hist1 and LAD columns)
    #clusters = list of clusters (each cluster is a list of NP names)
    results = []

    for cluster_id, cluster in enumerate(clusters, start=1):

        for np_name in cluster:

            total_windows_in_np = 0
            hist1_count = 0
            lad_count = 0

            for i in range(len(df)):

                #if this window is present in the NP
                if df[np_name].iat[i] == 1:

                    total_windows_in_np += 1

                    #check Hist1
                    if feat["Hist1"].iat[i] == 1:
                        hist1_count += 1

                    #check LAD
                    if feat["LAD"].iat[i] == 1:
                        lad_count += 1

            #avoid division by zero
            if total_windows_in_np > 0:
                hist1_pct = 100 * (hist1_count / total_windows_in_np)
                lad_pct = 100 * (lad_count / total_windows_in_np)
            else:
                hist1_pct = 0
                lad_pct = 0

            #puts results into a df
            results.append({ 
                "cluster": cluster_id,
                "NP": np_name,
                "Hist1_pct": hist1_pct,
                "LAD_pct": lad_pct
            })

    return pd.DataFrame(results)

def plot_feature_boxplots(results_df, value_col, ylabel, title, outpath):

    data = [
        results_df[results_df["cluster"] == 1][value_col],
        results_df[results_df["cluster"] == 2][value_col],
        results_df[results_df["cluster"] == 3][value_col]
    ]

    plt.figure(figsize=(8,5))
    plt.boxplot(data, labels=["Cluster 1", "Cluster 2", "Cluster 3"])
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close()











######## ACTIVITY 2 FUNCTIONS ###########



######## ACTIVITY 4 FUNCTIONS ###########

#computes the percentage of each feature for each NP in each cluster, returns a df with cluster, NP, feature, and percent columns
def compute_all_feat_stats(df, feat, clusters, desired_feats: list[str]):
    
    results = []
    
    for index, cluster in enumerate(clusters, start=1):
        for np_name in cluster:
            for feature in desired_feats:
                feature_count = 0
                total_windows = 0
                
                for i in range(len(df)):
                    # total_windows += 1
                    if df.iloc[i][np_name] == 1:
                        total_windows += 1
                    if df.iloc[i][np_name] == 1 and feat.iloc[i][feature] == 1:
                        feature_count += 1
                if total_windows > 0:
                    feat_percent = (feature_count / total_windows) * 100
                else:
                    feat_percent = 0

                    #puts results into a df
                results.append({ 
                    "cluster": index,
                    "NP": np_name,
                    "feature": feature,
                    "percent": feat_percent #percent is the percent of windows present in that np have features matching the feature in the feature table

                })

    return pd.DataFrame(results)


#creates radar chart for each cluster showing the average percent for each feature in that cluster
#useful website: https://medium.com/top-python-libraries/how-to-draw-radar-chart-in-python-3f56e52bcb17
def plot_radar(cluster_id, row):
    features = row.index.tolist()
    values = row.values.tolist()

    angles = np.linspace(0, 2*np.pi, len(features), endpoint=False).tolist()

    #close circle
    values += values[:1]
    angles += angles[:1]


    fig, ax = plt.subplots(figsize=(8,8), subplot_kw=dict(polar=True))
    ax.plot(angles, values)
    ax.fill(angles, values, alpha=0.15)

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(features)
    y_max = max(values)
    ax.set_ylim(0, y_max)

    ax.set_title(f"Cluster {cluster_id} Feature Profile")
    plt.tight_layout()

    os.makedirs("./radar_charts", exist_ok=True)
    plt.savefig(f"./radar_charts/cluster_{cluster_id}_radar.png")
    plt.close()
                
        


######## ACTIVITY 4 FUNCTIONS ###########


############## --- co-segregation functions --- ###############

#calculates the detection frequency for a given window (the fraction of NPs that have that window)
def calc_detection_freq(window):
    freq = window.mean()
    return freq

#calculates the co-segregation frequency for two windows (the fraction of NPs that have both windows)
def calc_cosegregation(window_a, window_b, np_list):
    return (((window_a == 1) & (window_b == 1)).sum() / len(np_list))
    # return ((window_a == 1) & (window_b == 1)).mean() #fraction of NPs that have both windows

#calculates the difference between ovserved co-segregation and expected co-segregation
def linkage(freq_a, freq_b, freq_ab):
    return freq_ab - (freq_a * freq_b) #the observed co-segregation minus the expected co-segregation 

#calculates the normalized linkage for two windows
def normalized_linkage(freq_a, freq_b, freq_ab):
    D = linkage(freq_a, freq_b, freq_ab) #the difference between observed and expected co-segregation

    #dmax is the maximum possible value of D based on the freqs of the two windows
    if D == 0:
        return 0
    elif D > 0:
        dmax = min(freq_a * (1 - freq_b), freq_b * (1 - freq_a))
    elif D < 0:
        dmax = min(freq_a * freq_b, (1 - freq_a) * (1 - freq_b))
    if dmax == 0:
        return 0

    return D / dmax

#builds the normalized linkage matrix
def compute_normalized_linkage_matrix(df, np_cols):
    n = len(df) #number of rows/windows in the df
    
    #creates the labels for the windows in the format "chr:start-stop"
    window_labels = (
        df["chrom"] + ":" +
        df["start"].astype(str) + "-" +
        df["stop"].astype(str)
    ).tolist()

    #initializes a dataframe of zeros of size num_windows x num_windows to store the normalized linkage values
    norm_linkage_matrix = pd.DataFrame(
        np.zeros((n, n)),
            index=window_labels, 
            columns=window_labels
    )


    frequencies = [] #all of the detection frequencies for each window

    #calculates the freqquencies and populates the frequencies list
    for i in range(n):
        window = df.iloc[i][np_cols]
        frequencies.append(calc_detection_freq(window))

    #loops through all pairs of windows to calculate the normalized linkage and populate the matrix
    x = 0
    for i in range(n):
        for j in range(i, n):
            window_a = df.iloc[i][np_cols] #grabs one of the windows as a series of 0s and 1s for each NP
            window_b = df.iloc[j][np_cols] 

            freq_a = frequencies[i]
            freq_b = frequencies[j]
            freq_ab = calc_cosegregation(window_a, window_b, np_cols) #gets co-seg for the pair of windows


            norm_linkage = normalized_linkage(freq_a, freq_b, freq_ab) #calculates the normalized linkage for the pair of windows
            norm_linkage_matrix.iat[i, j] = norm_linkage #fills in the upper triangle and diagonal
            norm_linkage_matrix.iat[j, i] = norm_linkage #fills in the lower triangle with same values

    return norm_linkage_matrix 


############## --- co-segregation functions --- ###############



############## --- network functions --- ###############

#makes the adj matrix, where two windows are connected if their normalized linkage is above the 75th percentile (Q3) of all normalized linkage values
def create_adj_matrix(df, Q3):
    n = len(df)

    adj_matrix = pd.DataFrame(
        np.zeros((n, n)),
        index=df.index, 
        columns=df.index
    )

    for i in range(n):
        for j in range(i+1,n):
            if df.iloc[i,j] > Q3:
                adj_matrix.iloc[i,j] = 1
                adj_matrix.iloc[j,i] = 1
    
    return adj_matrix

#explain in just a small sentence
#returms a series of degree centrality values for each node in the graph, calculated as the number of connections divided by the maximum possible connections (n-1)
#tldr: returns an array of degree centrality values for each node in the graph, where degree centrality is the number of connections a node has divided by the maximum possible connections (n-1)
def degree_centrality(adj_matrix):
    n = len(adj_matrix)
    return (adj_matrix.sum(axis=1) / (n - 1)).sort_values()


def print_degree_centrality_stats(dc):
    print("Average Degree Centrality:", dc.mean())
    print("Min Degree Centrality:", dc.min())
    print("Max Degree Centrality:", dc.max())
    print(dc)


def network_graph(adj_matrix, out_path):
    import pandas as pd
    import networkx as nx
    import matplotlib.pyplot as plt

    G = nx.from_pandas_adjacency(adj_matrix)
    G.name = "Adjacency Graph"
    print(G)
    plt.figure(figsize=(12, 12))
    pos = nx.spring_layout(G, seed=42, k=0.5)
    nx.draw_networkx(G, pos=pos, with_labels=False, font_size=6, node_size=50, width=0.5)

    plt.title("Network Graph of Genomic Interactions")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()

############## --- network functions --- ###############





############## --- community detection functions --- ###############

# def single_degree_centrality(adj_matrix, node):
    

#grabs the top five nodes with the highest degree centrality values and returns them as a series with the node names as the index and the degree centrality values
def top_five_degree_centrality(dc):
    return dc.nlargest(5)

#clusters all nodes into one of the five clusters with the highest degree centrality
def cluster_by_top_degree_centrality(adj_matrix, linkage_matrix, top_nodes):
    n = len(adj_matrix)
    clusters = {node: [] for node in top_nodes} #initializes a dictionary with the top nodes as keys and empty lists as values to store the clusters

    for node in adj_matrix.index:
        if node in top_nodes: #matches the top nodes to their own cluster
            clusters[node].append(node)
        else:
            max_centrality = -1 #starts with a max centrality of -1 to ensure any real centrality will be higher
            assigned_cluster = None

            for top_node in top_nodes:
                if adj_matrix.loc[node, top_node] == 1: #if there is a connection
                    link_strength = linkage_matrix.loc[node, top_node] #gets the normalized linkage value for the connection
                    if link_strength > max_centrality: #if this connection has a higher centrality than the previous max, update the assigned cluster
                        max_centrality = link_strength
                        assigned_cluster = top_node

            if assigned_cluster is not None:
                clusters[assigned_cluster].append(node)

    return clusters


#The size of the community (number of nodes).
# • The percentage of nodes in the community that contain a hist1 gene.
# • The percentage of nodes in the community that contain a LAD.
# • The list of nodes that are in the community
def print_community_information(clusters, feat):
    for top_node, cluster_nodes in clusters.items():
        size = len(cluster_nodes)
        hist1_count = 0
        lad_count = 0

        for node in cluster_nodes:
            if feat.loc[feat["name"] == node, "Hist1"].values[0] == 1:
                hist1_count += 1
            if feat.loc[feat["name"] == node, "LAD"].values[0] == 1:
                lad_count += 1

        hist1_pct = (hist1_count / size) * 100 if size > 0 else 0
        lad_pct = (lad_count / size) * 100 if size > 0 else 0

        print(f"Community centered on {top_node}:")
        print(f"Size: {size}")
        print(f"Percentage with Hist1: {hist1_pct:.2f}%")
        print(f"Percentage with LAD: {lad_pct:.2f}%")
        print(f"Nodes in community: {cluster_nodes}\n")











############## --- community detection functions --- ###############


