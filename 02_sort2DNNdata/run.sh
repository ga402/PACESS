#!/bin/sh                                                                                                                                                    

# The purpose of this script is to generate a dataframe from 
# the output of 2D object-detection neural network which is 
# ordered and sorted appropriately for creating 3d predictions (ordered by size/intensity)
# e.g

#data_path = "myDataFrame.csv"
#image_path = "bonemarrow.npy"
#output_path = "output.csv"
#channgel_order = [1,2]
#color_order = ["dTomato", "gfp"]
#selected_colors =  ["dTomato"] (e.g. you only want to select this channel to output from as maybe GFP represents a structural elements)

# Variables:
INFILE=input.csv # the 

OUTFILE=output.csv
IMAGE="myimage.npy" # saved as numpy file

julia getIntensityData.jl --data $INFILE \
                        --image $IMAGE \
                        --output $OUTFILE \
                        --channel_order 1 2 3 4 \
                        --color_order "dapi" "gfp" "bone" "dTomato" --select_colors "gfp" "dTomato" "all"

