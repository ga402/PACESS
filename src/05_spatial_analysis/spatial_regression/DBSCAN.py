# This code is written for readability rather than efficiency.
# 
# - It is aimed to give a set of logical steps explaining how to run the analysis
# in the manner we ran it for this study. 
# 

# 1. load up required packages and read the data

import pandas as pd
import numpy as np
from scipy.spatial import cKDTree

# Load data
df = pd.read_csv("your_file.csv")  # Replace with your actual filename
df['id'] = np.arange(1, len(df) + 1)
df = df.iloc[:, [6] + list(range(6))]  
df['sumTF'] = 0

# 2. Neighbor aggregation function
def colse_function(df, thirdqu):
    df = df.copy()
    for i in range(len(df)):
        a = df.iloc[i]

        up_bottom = df[(df['x'] == a['x']) & (df['y'] == a['y']) & (df['z'].isin([a['z'] + 45, a['z'] - 45]))]
        left_right = df[(df['x'] == a['x']) & (df['z'] == a['z']) & (df['y'].isin([a['y'] + 45, a['y'] - 45]))]
        front_back = df[(df['y'] == a['y']) & (df['z'] == a['z']) & (df['x'].isin([a['x'] + 45, a['x'] - 45]))]

        b = pd.concat([pd.DataFrame([a]), up_bottom, left_right, front_back])

        if b['AML'].mean() >= thirdqu:
            df.loc[df['id'].isin(b['id']), 'sumTF'] += 1

    return df

## Apply function
newdf = colse_function(df, 13)
newdf['cluster_num'] = 0

# 3. Clustering function
def onecluster_function(k):
    global newdf
    cluster = newdf[(newdf['sumTF'] > 0) & (newdf['cluster_num'] == 0)]
    if cluster.empty:
        return newdf

    aid = [cluster.iloc[0]['id']]
    len_aid0, len_aid1 = 0, 1

    while len_aid1 > len_aid0:
        len_aid0 = len(aid)
        a = cluster[cluster['id'].isin(aid)]
        left_a = cluster[~cluster['id'].isin(aid)]

        if left_a.empty:
            break

        tree = cKDTree(left_a[['x', 'y', 'z']].values)
        distances, indices = tree.query(a[['x', 'y', 'z']].values, k=min(len(left_a), len(a)))

        # Ensure arrays are 2D
        if distances.ndim == 1:
            distances = distances[:, np.newaxis]
            indices = indices[:, np.newaxis]

        match_idx = np.argwhere(np.isclose(distances, 45))
        new_ids = set()

        for i, j in match_idx:
            new_ids.add(left_a.iloc[indices[i, j]]['id'])

        aid = list(set(aid).union(new_ids))
        len_aid1 = len(aid)

    newdf.loc[newdf['id'].isin(aid), 'cluster_num'] = k
    return newdf

## Cluster all components
k = 1
while not newdf[(newdf['sumTF'] > 0) & (newdf['cluster_num'] == 0)].empty:
    newdf = onecluster_function(k)
    k += 1

# 4. Summary function
def sum_function(newdf):
    num_AML, num_T, num_MGK, num = [], [], [], []

    for i in range(1, newdf['cluster_num'].max() + 1):
        a = newdf[newdf['cluster_num'] == i]
        num.append(len(a))
        num_AML.append(a['AML'].sum())
        num_T.append(a['TC'].sum())
        num_MGK.append(a['MGK'].sum())

    sum_data = pd.DataFrame({
        'id': range(1, newdf['cluster_num'].max() + 1),
        'num': num,
        'num_AML': num_AML,
        'num_MGK': num_MGK,
        'num_T': num_T
    })

    sum_data['cluster'] = sum_data['num_AML'].rank(ascending=False, method='first').astype(int)
    sum_data['per_AML'] = sum_data['num_AML'] / newdf['AML'].sum()

    return sum_data

## Apply summary
sum_data = sum_function(newdf)

# 5. Map cluster ranks
newdf['cluster'] = 0
for i in range(1, newdf['cluster_num'].max() + 1):
    newdf.loc[newdf['cluster_num'] == i, 'cluster'] = sum_data.loc[sum_data['id'] == i, 'cluster'].values[0]

# Show cluster distribution
print(newdf['cluster'].value_counts())

# 6. Export final data
finaldata = newdf.iloc[:, [1, 2, 3, 4, 5, 6, 9]]  
finaldata.to_csv("DBSCAN.csv", index=False)













# 2. Convert loc string to coordinates

coords = new['loc'].str.split(',', expand=True).astype(float)
coords.columns = ['x', 'y', 'z']
DD = coords.copy()


# 3. DBSCAN clustering

db = DBSCAN(eps=45, min_samples=(13 + 1) * 7)
clusters = db.fit_predict(DD)

DD['cluster'] = clusters
DD.to_csv("AML_DBSCAN1.csv", index=False)


# 4. Summarise location by cluster count

df = pd.read_csv("AML_DBSCAN1.csv")
df['loc'] = df[['x', 'y', 'z']].astype(str).agg(','.join, axis=1)

l = []
un = []
cluster = []

for loc_val in df['loc'].unique():
    subset = df[df['loc'] == loc_val]
    l.append(len(subset))
    un.append(subset['cluster'].nunique())
    cluster.append(subset['cluster'].iloc[0])

loc_df = pd.DataFrame(df['loc'].unique(), columns=['loc'])
DD2 = loc_df['loc'].str.split(',', expand=True).astype(float)
DD2.columns = ['x', 'y', 'z']
DD2['AML'] = l
DD2['cluster'] = cluster



# 5. Count AML per cluster

AMLcluster = DD2[DD2['cluster'] != -1].groupby('cluster')['AML'].sum()
DD2['AML'] = DD2['AML'] - 1



# 6. Merge with Tc data

df1 = pd.read_csv("CD8final.csv")
dff = pd.merge(df1, DD2, on=["x", "y", "z", "AML"], how="left")
dff.to_csv("AML_DBSCAN2.csv", index=False)



# 7. Distance and cluster merge functions


def distfunction(a, b):
    return cdist(a[['x', 'y', 'z']], b[['x', 'y', 'z']])

def newclusterfunction(DD, grid):
    cluster_ids = sorted(set(DD['cluster']) - {-1})
    merge_matrix = np.zeros((len(cluster_ids), len(cluster_ids)), dtype=int)
    id_to_idx = {cid: i for i, cid in enumerate(cluster_ids)}

    for i, cid1 in enumerate(cluster_ids):
        a = DD[DD['cluster'] == cid1]
        x_range = (a['x'].min() - grid, a['x'].max() + grid)
        y_range = (a['y'].min() - grid, a['y'].max() + grid)
        z_range = (a['z'].min() - grid, a['z'].max() + grid)

        DDrange = DD[
            (DD['x'].between(*x_range)) &
            (DD['y'].between(*y_range)) &
            (DD['z'].between(*z_range))
        ]

        nearby_clusters = set(DDrange['cluster']) - {cid1, -1}
        for cid2 in nearby_clusters:
            b = DD[DD['cluster'] == cid2]
            if (distfunction(a, b) == grid).any():
                i1, i2 = id_to_idx[cid1], id_to_idx[cid2]
                merge_matrix[i1, i2] = merge_matrix[i2, i1] = 1

    DD['newcluster'] = DD['cluster']
    visited = set()

    for i, cid1 in enumerate(cluster_ids):
        if cid1 in visited:
            continue
        connected = {cid1}
        queue = [cid1]
        while queue:
            current = queue.pop()
            idx_current = id_to_idx[current]
            for j, val in enumerate(merge_matrix[idx_current]):
                if val == 1:
                    cid2 = cluster_ids[j]
                    if cid2 not in connected:
                        connected.add(cid2)
                        queue.append(cid2)
        visited.update(connected)
        min_id = min(connected)
        DD.loc[DD['cluster'].isin(connected), 'newcluster'] = min_id

    cluster_map = {cid: idx for idx, cid in enumerate(sorted(DD['newcluster'].unique()))}
    DD['finalcluster'] = DD['newcluster'].map(cluster_map)
    return DD

# 8. Run merge and export...
DD = pd.read_csv("AML_DBSCAN2.csv")
newDD = newclusterfunction(DD, grid=45)
newDD = newDD.iloc[:, list(range(6)) + [-1]]
newDD.to_csv("AML_DBSCAN3.csv", index=False)


# 9. Final relabel and export


df = pd.read_csv("AML_DBSCAN3.csv")
cluster_sizes = df.groupby('finalcluster')['AML'].sum().sort_values(ascending=False)

order_map = {cid: i+1 for i, cid in enumerate(cluster_sizes.index)}
df['cluster'] = df['finalcluster'].map(order_map)

dfz = df[['x', 'y', 'z', 'AML', 'MGK', 'TC', 'cluster']]
dfz.columns = ["x", "y", "z", "AML", "MGK", "TC", "cluster"]
dfz.to_csv("final.csv", index=False)
