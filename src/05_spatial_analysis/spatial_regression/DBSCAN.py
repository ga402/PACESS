# This code is written for readability rather than efficiency.
# 
# - It is aimed to give a set of logical steps explaining how to run the analysis
# in the manner we ran it for this study. 
# 


# 1. Read csv, expand AML counts into locations

df = pd.read_csv("AML_input.csv").iloc[:, :4]
df['AML'] = df['AML'] + 1

df['loc'] = df[['x', 'y', 'z']].astype(str).agg(','.join, axis=1)

loc = []
for i, row in df.iterrows():
    loc.extend([row['loc']] * row['AML'])

new = pd.DataFrame({'loc': loc})


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
