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






