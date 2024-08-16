
# Example 3D plot of clusters 
# Author: George Adams
# Date: August 2024


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
import geopandas as gpd
import duckdb
from matplotlib.colors import ListedColormap
# functions
from utils.functions import *
plt.rcParams['figure.dpi']=160




# import data...
path_name_gdf = "./boundaries.geojson"
gdf = gpd.read_file(path_name_gdf)

dt = pd.read_csv('leukaemia_sample1.csv')


# DATA PREPROCESSING/WRANGLING
# convert from pixels to um
dt.loc[:, ['x', 'y', 'z']] = dt.loc[:, ['x', 'y', 'z']] / [0.41, 0.41, 5]


# create dataset with AML in quartiles and groups of CTLs and MGKs
dg = duckdb.sql("""
    SELECT x, y, z, AML, ntile(4) OVER (ORDER BY AML) -1 as AML_GROUPS,
    CASE 
        WHEN
            (MGK = 0 AND TC >0) THEN 'TC'
        WHEN
            (MGK > 0 AND TC >0) THEN 'TC_MGK'
        WHEN 
            (MGK > 0 AND TC = 0) THEN 'MGK'
        ELSE
            'OTHER'
        END AS CELL_GROUPS,
    from dt
    ORDER BY AML_GROUPS ASC
    """).to_df()


# get quartile thresholds
quartile_thresholds = dg.groupby('AML_GROUPS').agg({'AML':[min, max, 'count']}).reset_index().iloc[:, 2].to_list()

# convert TC/TC_MGK/MGK into dummy variables
dg = pd.get_dummies(dg, columns=['CELL_GROUPS', ])


# -------------------------------------------------------


# create colour maps 
norm = mcolors.Normalize(vmin=0, vmax=1)
linline = np.linspace(0, 1, 256)

# yellow, blue, green  - cmaps from scratch!
yellow_colors = [(1, 0.76, 0, norm(value)) for value in linline]  # variable alpha value
blue_colors = [(0, 0, 1, norm(value)) for value in linline] 
green_colors = [(35/255, 202/255, 138/255, norm(value)) for value in linline]  

blue_cmap = mcolors.ListedColormap(blue_colors)
green_cmap = mcolors.ListedColormap(green_colors)
yellow_cmap = mcolors.ListedColormap(yellow_colors)


# for red... update the Red cmap 
red_cmap = plt.cm.Reds
my_cmap = red_cmap(np.arange(red_cmap.N))
my_cmap[:,-1] = np.sqrt(np.linspace(0, 1, red_cmap.N))
red_my_cmap = ListedColormap(my_cmap)


# -----------------------------------------------------

# get contours boundaries to ensure fit data...
gdf = gdf.assign(ymax = lambda x: x[['geometry']].bounds.maxy)
gdf = gdf.assign(xmax = lambda x: x[['geometry']].bounds.maxx)
ygmax = gdf.ymax.max()
xgmax = gdf.xmax.max()


# ----------------------------------------------------

# plots...

with plt.style.context(['science']):

    green_color_ = (35/255, 202/255, 138/255,1)
    blue_color_ = (0, 0, 1,1)
    
    zlevels =np.sort(dg.z.unique())
    
    fig = plt.figure(figsize=(8,6))
    ax0 = fig.add_subplot(1, 1, 1, projection='3d')

    for n,z in enumerate(zlevels):
        dfz = dg.query(f'z=={z}')
        
        # plot the AML regions 
        arr = getBlurArray(dfz, z, 'AML_GROUPS', size=2)
        Ny, Nx = arr.shape
        XX, YY, ZZ = np.meshgrid(np.arange(Nx), np.arange(Ny), np.arange(z_max))
        arr_alpha = arr.copy()
        cntr1 = ax0.contourf(
            XX[:, :, 10], YY[:, :, 10], arr_alpha,
            zdir='z', offset=+z*5,cmap=red_my_cmap, zorder=0,levels=[0.5, 1.5, 2.5, 3.5, 4.5]
            )
        # add a color bar
        if z == 0:
            cbar = plt.colorbar(cntr1, 
                                orientation='horizontal',
                                ticks=[0.5, 1.5, 2.5, 3.5, 4.5],
                                format=mticker.FixedFormatter([0] + quartile_thresholds),
                                pad=0.0, shrink=0.6)
            cbar.set_label('AML quartiles ranges', rotation=0, labelpad=10)
        
        
        # plot MGK/CTL groups
        for gr,cmt in zip(['CELL_GROUPS_MGK', 'CELL_GROUPS_TC'], [blue_cmap, green_cmap]):
            arr = getBlurArray(dfz, z, gr, size=3)
            Ny, Nx = arr.shape
            XX, YY, ZZ = np.meshgrid(np.arange(Nx), np.arange(Ny), np.arange(z_max))
            arr_alpha = arr.copy()
            _ = ax0.contourf(
                XX[:, :, 10], YY[:, :, 10], arr_alpha
                zdir='z', offset=+z*5,cmap=cmt, zorder=0
            )
        
        # plot bone contours
        for i in range(list(gdf.query(f"z_level=={z+4}").geometry.iloc[0]).__len__()):
            xx, yy = zip(*list(gdf.query(f"z_level=={z+4}").geometry.iloc[0])[i].exterior.coords)
            x = np.array(xx) * np.array(Nx*1.2 /xgmax) -10
            y = np.array(yy) * np.array(Ny /ygmax)
            for ax in [ax0]:
                ax.plot(x, y, z*5, c='k', zorder = 5+ z)
  
    
    for ax in [ax0]:
        labs = ax.get_xticks()
        ax.set_xticklabels([roundup(l*45) if n&1 else '' for n,l in enumerate(labs)])
        labs = ax.get_yticks()
        ax.set_yticklabels([roundup(l*45) if n&1 else '' for n,l in enumerate(labs)])
        labs = ax.get_zticks()
        ax.set_zticklabels([int(round(float(l),0)) for l in labs])
        ax.view_init(elev=22, azim=-22)
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False

        # Now set color to white (or whatever is "invisible")
        ax.xaxis.pane.set_edgecolor('w')
        ax.yaxis.pane.set_edgecolor('w')
        ax.zaxis.pane.set_edgecolor('w')
        green_patch = mpatches.Patch(color=green_color_, label='CTL locations')
        blue_patch = mpatches.Patch(color=blue_color_, label='MGK locations')
        ax.legend(handles=[green_patch, blue_patch])


    plt.tight_layout()
    # set the spacing between subplots
    plt.subplots_adjust(left=0.2,
                        bottom=0, 
                        right=0.8, 
                        top=1, 
                        wspace=0, 
                        hspace=0)
    

    plt.savefig('../../figure4G.png')
















