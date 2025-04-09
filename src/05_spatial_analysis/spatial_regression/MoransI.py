

# This code is written for readability rather than efficiency.
# 
# - It is aimed to give a set of logical steps explaining how to run the analysis
# in the manner we ran it for this study. 
# 



import pandas as pd
import numpy as np
from distance3d import distance3d  # Cython version of rcpp_distance3d

# 1. Read data and rename columns
df = pd.read_csv("your_file.csv", header=None)
df.columns = ["x", "y", "z", "AML", "MGK", "Tcell"]

# 2. Compute 3D Euclidean distance matrix using your Cython function
coords = df[['x', 'y', 'z']].values
dmat = distance3d(coords)

# 3. Moran’s I Function (binary kernel at distance == 45)
def morans_i(a, dmat):
    dmat1 = (dmat == 45).astype(int)

    m = np.zeros_like(a, dtype=float)
    n = np.zeros(len(a), dtype=float)
    n1 = np.zeros(len(a), dtype=float)

    for i in range(len(a)):
        m_vec = np.where(
            dmat1[i] == 1,
            (a[i] - np.mean(a)) * (a - np.mean(a)),
            0
        )
        n[i] = np.sum(m_vec)
        n1[i] = (a[i] - np.mean(a)) ** 2

    I = (len(a) * np.sum(n)) / (np.sum(dmat1) * np.sum(n1))
    return I

# 4. Apply for each cell type
I_tcell = morans_i(df["Tcell"].values, dmat)
I_aml = morans_i(df["AML"].values, dmat)
I_mgk = morans_i(df["MGK"].values, dmat)

# 5. Collect and save
I = pd.Series([I_tcell, I_aml, I_mgk], index=["Tcell", "AML", "MGK"])
I.to_csv("MoransI.csv", header=False)
