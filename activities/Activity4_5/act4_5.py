import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import random

path = "../data/GSE64881_segmentation_at_30000bp.passqc.multibam.txt"
df2 = pd.read_csv(path, sep='\t') #hist1 region


# range_to_keep = [21700000, 24100000]

#Creates a deep copy in a dataframe for chrom 13 between the desired coordinates
df2 = df2[(df2['chrom'] == 'chr13') & (df2['start'] >= 21690000) & (df2['stop'] <= 24120000)].copy()


#goes through and eliminates the columns without any windows present
for col in df2.columns[3:]:
    if df2[col].sum() == 0:
        df2 = df2.drop(columns=col)
        

                 
np_cols = df2.columns[3:] #useable np columns, ignoring the first three cols

#the matrix to be filled
jaccard_matrix = pd.DataFrame(
    np.eye(len(np_cols)), #ones along the main diagonal
    index = np_cols, #assigns row labels
    columns = np_cols #assigns col labels
)


#populates the jaccard index matrix
for a in range(len(np_cols)): 
    for b in range(a, len(np_cols)): #start at a until end
        m_00 = m_01 = m_10 = m_11 = 0
        
        colA = df2[np_cols[a]] 
        colB = df2[np_cols[b]]
        
        for x in range(len(df2)):
            if colA.iat[x] == 1 and colB.iat[x] == 1:
                m_11 += 1
            elif colA.iat[x] == 1 and colB.iat[x] == 0:
                m_10 += 1
            elif colA.iat[x] == 0 and colB.iat[x] == 1:
                m_01 += 1
        
        if ((m_01 + m_10 + m_11) > 0): #makes sure its not zero to avoid division by zero
            jaccard = m_11 / (m_01 + m_10 + m_11)
        else:
            jaccard = 0
            
        jaccard_matrix.iat[a,b] = jaccard
        jaccard_matrix.iat[b,a] = jaccard
        
        
dist_matrix = 1 - jaccard_matrix #makes the dist_matrix


########### Heatmap output for jaccard index and jaccard distance matracies ########


# plt.figure(figsize=(12, 10))

# sns.heatmap(
#     jaccard_matrix,
#     cmap="viridis",
#     square=True,
#     vmin=0,
#     vmax=1,
#     xticklabels=False,
#     yticklabels=False,
#     cbar_kws={"label": "Jaccard index"}
# )

# plt.title("Jaccard Similarity Heatmap (Hist1 Region)")
# plt.tight_layout()

# plt.savefig(
#     "jaccard_index_heatmap.png",
#     dpi=300,
#     bbox_inches="tight"
# )

# plt.close()

# plt.figure(figsize=(12, 10))

# sns.heatmap(
#     dist_matrix,
#     cmap="viridis",
#     square=True,
#     vmin=0,
#     vmax=1,
#     xticklabels=False,
#     yticklabels=False,
#     cbar_kws={"label": "Jaccard distance"}
# )

# plt.title("Jaccard Distance Heatmap (Hist1 Region)")
# plt.tight_layout()

# plt.savefig(
#     "jaccard_distance_heatmap.png",
#     dpi=300,
#     bbox_inches="tight"
# )

# plt.close()


########### Heatmap output for jaccard index and jaccard distance matracies #########


#------------ Clustering Part 1.  ------------#

#dataframe of size np_cols x np_cols
norm_jaccard_matrix = pd.DataFrame(
    np.eye(len(np_cols)),
    index=np_cols, #assigns row labels
    columns=np_cols #assigns col labels
)

num_ones_a = num_ones_b = m_11 = 0


#populates the normalized jaccard index matrix
for a in range(len(np_cols)): 
    for b in range(a, len(np_cols)): #start at 'a' until end to prevent double counting
        num_ones_a = num_ones_b = m_11 = 0
        
        colA = df2[np_cols[a]] 
        colB = df2[np_cols[b]]
        
        for x in range(len(df2)):
            if colA.iat[x] == 1 and colB.iat[x] == 1:
                m_11 += 1
            if colA.iat[x] == 1:
                num_ones_a += 1
            if colB.iat[x] == 1:
                num_ones_b += 1
                
        denom = min(num_ones_a,num_ones_b)
        
        if ((denom) > 0): #makes sure its not zero to avoid division by zero
            jaccard_sim = m_11 / denom
        else:
            jaccard_sim = 0
            
        norm_jaccard_matrix.iat[a,b] = jaccard_sim
        norm_jaccard_matrix.iat[b,a] = jaccard_sim
        
norm_distance_matrix = 1 - norm_jaccard_matrix #makes the normalized distance matrix



three_points = np.random.choice(np_cols, size=3, replace=False) 

print(three_points) #the three random points


dist_to_centers = norm_distance_matrix[three_points] #163x3 dataframe



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
            


clusters = assign_to_k_clusters(dist_to_centers, three_points)
print("Cluster 1:",len(clusters[0]))
print("Cluster 2:",len(clusters[1]))
print("Cluster 3:",len(clusters[2]))

#------------ Clustering Part 1.  ------------#





######## Normalized heat maps ##########

# plt.figure(figsize=(12, 10))

# sns.heatmap(
#     norm_jaccard_matrix,
#     cmap="viridis",
#     square=True,
#     vmin=0,
#     vmax=1,
#     xticklabels=False,
#     yticklabels=False,
#     cbar_kws={"label": "Normalized Jaccard index"}
# )

# plt.title("Normalized Jaccard Similarity Heatmap (Hist1 Region)")
# plt.tight_layout()

# plt.savefig(
#     "norm_jaccard_index_heatmap.png",
#     dpi=300,
#     bbox_inches="tight"
# )

# plt.close()





# plt.figure(figsize=(12, 10))

# sns.heatmap(
#     norm_distance_matrix,
#     cmap="viridis",
#     square=True,
#     vmin=0,
#     vmax=1,
#     xticklabels=False,
#     yticklabels=False,
#     cbar_kws={"label": "Normalized Distance Matrix"}
# )

# plt.title("Normalized Distance Matrix Heatmap (Hist1 Region)")
# plt.tight_layout()

# plt.savefig(
#     "norm_distance_heatmap.png",
#     dpi=300,
#     bbox_inches="tight"
# )

# plt.close()

######## Normalized heat maps ##########




#############------ Cluster Part 2. -------###############


#finds and returns the np in a cluster that has the minimum disssimilarity
def find_center(cluster: list[str], norm_dist_matrix: pd.DataFrame) -> str:

    #best_np = cluster[0]
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
    
    
    
#calculates the variance in a cluster  
# def cluster_var(cluster: list[str], center, norm_dist_matrix: pd.DataFrame):
#     if len(cluster) == 0:
#         return 0
    
#     distances = []
    
#     for np_name in cluster:
#         d = norm_dist_matrix.loc[np_name, center]
#         distances.append(d)
    
#     avg = sum(distances) / len(distances)
    
    
#     var = 0
#     for d in distances:
#         var += (d - avg) ** 2
    
#     var = var / len(distances)

#     return var


#sums up the distances between all points in a cluster to the center, and sums for all clusters
def within_cluster_var(clusters: list[list[str]], centers: list[str], norm_dist_matrix: pd.DataFrame):
    total = 0.0 #the total distances from centers, within cluster variation
    
    for i, cluster in enumerate(clusters): #goes through all of the clusters
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
        for cluster in clusters:
            new_centers.append(find_center(cluster, norm_dist_matrix))
        
        new_clusters = assign_to_k_clusters(norm_dist_matrix[new_centers], new_centers) #reclusters based on new centers
        wc_var = within_cluster_var(new_clusters, new_centers, norm_dist_matrix) #gets the within cluster var for the new clusters
        
        #output
        print("\nIteration:", it+1)
        print("Variation:",wc_var)
        print(new_centers, '\n')

        if new_centers == centers:
            break
        
        centers = new_centers
        clusters = new_clusters
    
    return new_centers, new_clusters #returns the updated centers and clusters         
 


var_for_original = within_cluster_var(clusters, three_points, norm_distance_matrix)
print("Variation for original points:",var_for_original) #within cluster variation for the original random points



new_centers, clusters = cluster_medoids(norm_distance_matrix, three_points)


wc_var = within_cluster_var(clusters, new_centers, norm_distance_matrix)



print("Final Total Within Cluster Variation:", wc_var)

print("Final Centers:", new_centers)

print("Cluster 1:",len(clusters[0]))
print("Cluster 2:",len(clusters[1]))
print("Cluster 3:",len(clusters[2]))

#############------ Cluster Part 2. -------###############
