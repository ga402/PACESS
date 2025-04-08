#!/bin/sh                                                                                      

# This script is to generate the 3D predictions
# The output of this script is a .csv file for **EACH** different cell type for which you want a prediction for
# e.g. OUTPUTc2.csv, OUTPUTc4.csv, etc...
# Each of these files is the 3D predictions for this cell type.ss


INFILE=input.csv
OUTFILE=output.csv
UTILS_FILE=getDistancesUtils.jl

python main.py --file $INFILE --output $OUTFILE --pixelsXYZ 0.87 0.87 5 --cells 0 2 4 --intensity_match c2 c4 c0




