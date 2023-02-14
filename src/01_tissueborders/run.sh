#!/bin/zsh

# The purpose of this script is generate a 3D contour map (boundary map) of the tissue 
# in question. This will be used in later scripts to exclude (?anomolous) cells which
# are identified outside the boundary of the tissue and to generate the final grid 
# for the spatial dataframe.


$file=INPUT.tiff
$output=OUTPUT.geojson # geopandas dataframe of tissue contours 
python createBoundaries.py --file $file --output $output
