
# This code is written for readability rather than efficiency.
# 
# - It is aimed to give a set of logical steps explaining how to run the analysis
# in the manner we ran it for this study. 
# 



import pandas as pd
import numpy as np
from scipy.spatial.distance import cdist

# Load data
df = pd.read_csv("your_file.csv")
df = df.iloc[:, 1:]  # Drop the first column

# Distance function (3D Euclidean matrix between two sets of points)
def distfunction(a, b):
    coords_a = a[['x', 'y', 'z']].values
    coords_b = b[['x', 'y', 'z']].values
    return cdist(coords_a, coords_b)

# Main function to calculate TC and MGK differences
def difffunction(df):
    cluster_ids = sorted(df['cluster'].unique())
    cluster_ids = [cid for cid in cluster_ids if cid > 0]  # exclude background

    TCdiff = []
    MGKdiff = []

    for i in cluster_ids:
        a = df[df['cluster'] == i]

        # Define neighborhood cube around cluster i
        minx, maxx = a['x'].min() - 90, a['x'].max() + 90
        miny, maxy = a['y'].min() - 90, a['y'].max() + 90
        minz, maxz = a['z'].min() - 90, a['z'].max() + 90

        b = df[
            (df['x'] > minx) & (df['x'] < maxx) &
            (df['y'] > miny) & (df['y'] < maxy) &
            (df['z'] > minz) & (df['z'] < maxz)
        ]

        b = b[b['cluster'] != i]

        if b.empty:
            TCdiff.append(np.nan)
            MGKdiff.append(np.nan)
            continue

        dist = distfunction(a, b)
        mask = (dist == 45)

        near_b = b.iloc[np.where(mask.any(axis=0))[0]]

        if near_b.empty:
            TCdiff.append(np.nan)
            MGKdiff.append(np.nan)
            continue

        # Tag cluster 0 for "outside" and cluster >0 for "inside"
        near_b = near_b.copy()
        near_b['cluster'] = 0
        a['cluster'] = i  # already > 0

        final = pd.concat([a, near_b])

        tc0 = final.loc[final['cluster'] == 0, 'TC'].mean()
        tc1 = final.loc[final['cluster'] > 0, 'TC'].mean()
        mgk0 = final.loc[final['cluster'] == 0, 'MGK'].mean()
        mgk1 = final.loc[final['cluster'] > 0, 'MGK'].mean()

        TCdiff.append(tc0 - tc1)
        MGKdiff.append(mgk0 - mgk1)

    return pd.DataFrame({'TCdiff': TCdiff, 'MGKdiff': MGKdiff})

# Run and save
diff = difffunction(df)
diff.to_csv("diff.csv", index=False)
