import pandas as pd
import numpy as np







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




def load_and_filter(path):
    # loads the file and keeps just the hist1 region (chr13 range we used)
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
        cmap="viridis",
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

