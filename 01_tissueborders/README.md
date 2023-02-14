# Generate Image borders/extract contours

The purpose of the following script is to generate a geopandas dataframe of contours for each Z level of your image. 

| z_level | geometry |
|:-------:|:--------:|
| 1| [[0.1, 0.2, 0.1], ...] |
| 2 | [[0.1, 0.4,0.2],...] |
| ... | ...|
| 75 | [0.7, 0.2, 0.1], ...] |

To run this script, simply use the run.sh script and change the input and output parameters. 

The input file should be a binary tiff file in which your tissue has been roughly segmented to identify the borders you wish to extract. 
