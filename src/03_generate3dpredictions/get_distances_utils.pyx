# cython: boundscheck=False
import numpy as np
import pandas as pd
import re
from scipy.spatial import KDTree
from libc.math cimport fmax
from cython.parallel import prange

def getZlevel(path_name):
    return int(re.search(r"Z(\d+)", path_name).group(1))

def getXlevel(path_name):
    return int(re.search(r"X(\d+)", path_name).group(1))

def getYlevel(path_name):
    return int(re.search(r"Y(\d+)", path_name).group(1))

def addBegins(df):
    df['X_begin'] = df['image_name'].apply(getXlevel)
    df['Y_begin'] = df['image_name'].apply(getYlevel)
    df['Z_begin'] = df['image_name'].apply(getZlevel)
    return df

def sortMaxArea(df_m):
    df_m['maxArea'] = df_m['w'] * df_m['h']
    return df_m.sort_values('maxArea', ascending=False)

def sortMaxArea_Intensity(df, dict_cells, cell):
    df['maxArea'] = df['w'] * df['h']
    intensity_col = f"{dict_cells[cell]}Intensity"
    return df.sort_values([intensity_col, 'maxArea'], ascending=[False, False])

def calcMFI(img, df):
    assert img.ndim == 3
    cdef Py_ssize_t i, n = len(df)
    v = np.zeros(n, dtype=np.float64)

    for i in range(n):
        z = int(df.iloc[i]['z']) + 1
        y1 = int(df.iloc[i]['y'] + df.iloc[i]['Y_begin']) + 1
        y2 = int(df.iloc[i]['y'] + df.iloc[i]['Y_begin'] + df.iloc[i]['h'])
        x1 = int(df.iloc[i]['x'] + df.iloc[i]['X_begin']) + 1
        x2 = int(df.iloc[i]['x'] + df.iloc[i]['X_begin'] + df.iloc[i]['w'])
        region = img[z, y1:y2, x1:x2]
        v[i] = np.mean(region)
    return v

def createMatrix(df_m, x_c, y_c, z_c):
    mat = np.stack([
        df_m['x'].to_numpy() * x_c,
        df_m['y'].to_numpy() * y_c,
        df_m['z'].to_numpy() * z_c
    ])
    return mat.astype(np.float64)

def getDistancesBelowMaxDiam(ddf, x_c, y_c, z_c):
    mat = createMatrix(ddf, x_c, y_c, z_c).T
    tree = KDTree(mat)

    ddf['UpdateValue'] = 0
    taken_pos = set()

    for j in range(len(ddf)):
        r = max(ddf.iloc[j]['w'] * x_c, ddf.iloc[j]['h'] * y_c)
        if r > 250:
            continue
        idxs = tree.query_ball_point(mat[j], r)
        candidate_pos = set(idxs) - taken_pos
        if not candidate_pos:
            continue
        for idx in candidate_pos:
            ddf.at[idx, 'UpdateValue'] = ddf.iloc[j]['CellValueBase']
        taken_pos.update(candidate_pos)
    return ddf

def determineFinalCellPos(df_out):
    final_values = [
        row['UpdateValue'] if row['UpdateValue'] != 0 else row['CellValueBase']
        for _, row in df_out.iterrows()
    ]
    df_out['FinalCell'] = final_values
    return df_out

def identityEuclideanGroups(df, x_c, y_c, z_c, dict_cells, cell):
    df = sortMaxArea_Intensity(df, dict_cells, cell)
    df['CellValueBase'] = range(1, len(df)+1)
    df['UpdateValue'] = 0

    df = getDistancesBelowMaxDiam(df, x_c, y_c, z_c)
    df = determineFinalCellPos(df)
    return df
