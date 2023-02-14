#!/bin/zsh

# This script is to conbine all the 3D predictions into a 
# single spatial dataframe in which the 3d tissue space
# is divided into a set of n x n x n grids (size $n). In
# each grid is the COUNT/Number of each cell-type.

# to generate this grid - the contours of the 3d tissue are used (mygdpdata.geojson)

set -e

n=45 # the distance of the grid you want (um)

python generateRegularBoxCounts.py \
       --files #e.g. celltype1.csv celltype2.csv celltype3.csv... 
       --cells #e.g. leukaemia megakaryocyte cd8tcells \
       --xyz_correction # e.g. 0.414321 0.414321 5 \
       --edge_length $n \
       --gdf #e.g. mygpddata.geojson \
       --output #e.g.  myoutput$n_umGrid.csv \
       --save_data_out False # option
       --data_out # extra dataframe if above set to TRUE \


exit
