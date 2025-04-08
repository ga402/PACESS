#!/user/bin/env python3

import polars as pl
import numpy as np
from PIL import Image
import glob
from utils.image_loader import *
from utils.subsample_mean import *
import re
import logging
import argparse
from datetime import datetime
import sys


def get_program_parameters():
    epilogue = """
    get pixel intensities for the image
    """

    parser = argparse.ArgumentParser(
        description="sort and add labels",
        epilog=epilogue,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("-f", "--file", help="input .csv")
    parser.add_argument("-i", "--images", help= "input image folder")
    parser.add_argument("-o", "--output", help="output file")
    args = parser.parse_args()
    return args.file, args.images, args.output


def getZlevel(path_name:str) -> int:
    x = re.search(r"Z\d+", path_name).group(0)
    x = x.split('Z')[1]
    return int(x)


def getXlevel(path_name:str) -> int:
    x = re.search(r"X\d+", path_name).group(0)
    x = x.split('X')[1]
    return int(x)


def getYlevel(path_name:str) -> int:
    x = re.search(r"Y\d+", path_name).group(0)
    x = x.split('Y')[1]
    return int(x)






if __name__ == "__main__":
    logfile = datetime.now().strftime("logfile%Hhr%Mmin%d%m%Y.log")
    logging.basicConfig(filename=logfile, level=logging.DEBUG)


    input_path, image_path output_path = get_program_parameters()

    logging.info(f"{output_path=}")
    logging.info(f"{image_path=}")
    logging.info(f"{input_path=}")

    df = pl.read_csv(input_path)
    
    image_list = glob.glob(f'{image_path}/*.jpg')

    image_dict = load_images_to_dict(image_list)

	df = df.with_columns(
		pl.col("image_name").map_elements(getZlevel, return_dtype=pl.Int32).alias('Z_begin'),
		pl.col("image_name").map_elements(getXlevel, return_dtype=pl.Int32).alias('X_begin'),
		pl.col("image_name").map_elements(getYlevel, return_dtype=pl.Int32).alias('Y_begin')
		)


	df = df.with_columns(
		(pl.col('x') - pl.col('X_begin')).alias('x_pos'),
		(pl.col('y') - pl.col('Y_begin')).alias('y_pos')
		)



	df = df.with_columns(
		AML_MFI = vectorized_calc_mean_from_image_dict(image_dict, 
										df['image_name'].to_numpy(),
										df['x_pos'].to_numpy().astype(np.uint64),
										df['y_pos'].to_numpy().astype(np.uint64),
										df['w'].to_numpy().astype(np.uint64),
										df['h'].to_numpy().astype(np.uint64),
										int(0)
										),
		MGK_MFI = vectorized_calc_mean_from_image_dict(image_dict, 
										df['image_name'].to_numpy(),
										df['x_pos'].to_numpy().astype(np.uint64),
										df['y_pos'].to_numpy().astype(np.uint64),
										df['w'].to_numpy().astype(np.uint64),
										df['h'].to_numpy().astype(np.uint64),
										int(4)
										),
		TC_MFI = vectorized_calc_mean_from_image_dict(image_dict, 
									df['image_name'].to_numpy(),
									df['x_pos'].to_numpy().astype(np.uint64),
									df['y_pos'].to_numpy().astype(np.uint64),
									df['w'].to_numpy().astype(np.uint64),
									df['h'].to_numpy().astype(np.uint64),
									int(1)
									),

	)


    df.write_csv(f"{output_path}", separator=",")

    sys.exit()






# image_list = glob.glob('Images512x512_Z46/*.jpg')

# image_dict = load_images_to_dict(image_list)


# #gdf = gpd.read_file('./geodf_GEOM.geojson')

# df = pl.read_csv('SORTED_COMBINED_NN_OUTPUT_DATAFRAME_run5_noduplicates.csv')

# df = pl.read_csv('COMBINED_NN_OUTPUT_DATAFRAME_run5_noduplicates.csv')


# def getZlevel(path_name:str) -> int:
#     x = re.search(r"Z\d+", path_name).group(0)
#     x = x.split('Z')[1]
#     return int(x)


# def getXlevel(path_name:str) -> int:
#     x = re.search(r"X\d+", path_name).group(0)
#     x = x.split('X')[1]
#     return int(x)


# def getYlevel(path_name:str) -> int:
#     x = re.search(r"Y\d+", path_name).group(0)
#     x = x.split('Y')[1]
#     return int(x)



# df = df.with_columns(
# 	pl.col("image_name").map_elements(getZlevel, return_dtype=pl.Int32).alias('Z_begin'),
# 	pl.col("image_name").map_elements(getXlevel, return_dtype=pl.Int32).alias('X_begin'),
# 	pl.col("image_name").map_elements(getYlevel, return_dtype=pl.Int32).alias('Y_begin')
# 	)




# df = df.with_columns(
# 	(pl.col('x') - pl.col('X_begin')).alias('x_pos'),
# 	(pl.col('y') - pl.col('Y_begin')).alias('y_pos')
# 	)





# df = df.with_columns(
# 	AML_MFI = vectorized_calc_mean_from_image_dict(image_dict, 
# 									df['image_name'].to_numpy(),
# 									df['x_pos'].to_numpy().astype(np.uint64),
# 									df['y_pos'].to_numpy().astype(np.uint64),
# 									df['w'].to_numpy().astype(np.uint64),
# 									df['h'].to_numpy().astype(np.uint64),
# 									int(0)
# 									),
# 	MGK_MFI = vectorized_calc_mean_from_image_dict(image_dict, 
# 									df['image_name'].to_numpy(),
# 									df['x_pos'].to_numpy().astype(np.uint64),
# 									df['y_pos'].to_numpy().astype(np.uint64),
# 									df['w'].to_numpy().astype(np.uint64),
# 									df['h'].to_numpy().astype(np.uint64),
# 									int(4)
# 									),
# 	TC_MFI = vectorized_calc_mean_from_image_dict(image_dict, 
# 								df['image_name'].to_numpy(),
# 								df['x_pos'].to_numpy().astype(np.uint64),
# 								df['y_pos'].to_numpy().astype(np.uint64),
# 								df['w'].to_numpy().astype(np.uint64),
# 								df['h'].to_numpy().astype(np.uint64),
# 								int(1)
# 								),

# 	)

# df.head()



