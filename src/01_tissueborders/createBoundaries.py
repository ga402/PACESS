import numpy as np
import tifffile as tf

# import matplotlib.pyplot as plt
import cv2 as cv
import geopandas as gpd
from utils.contourtools import *
from utils.func import *
from utils.img2df import img2df as img2df


def get_program_parameters():
    epilogue = """
    create a dataframe of pixel sizes and z levels and image names
    """

    parser = argparse.ArgumentParser(
        description="create the dataframe",
        epilog=epilogue,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("-f", "--file", help="image path, control.tif")
    parser.add_argument(
        "-o", "--output", help="output name: eg. control_boundaries.geojson"
    )
    args = parser.parse_args()
    return args.file, args.output


@spapply3D
def main(im):
    """Main function for creating contours  - parallel function.

    Args:
        im ([type]): binary (padded) image (3D)

    Returns:
        list[numpy.ndarray]: contours list (3D)
    """
    contours, hierarchy = findContours(im)
    hierarchy = hierarchy.squeeze()
    filtered_contours, filtered_hierarchy = filterContoursSize(contours, hierarchy)
    filtered_contours2 = [
        c for c, h in zip(filtered_contours, filtered_hierarchy) if h[3] == -1
    ]
    contours_out = applySmoothContour(filtered_contours2, smooth_parameter=10)
    return contours_out


if __name__ == "__main__":
    logfile = createlog()
    logging.basicConfig(filename=logfile, level=logging.INFO)

    # image_path = "control.tif"
    # output_name = "control_boundaries.geojson"
    image_path, output_name = get_program_parameters()

    img = tf.imread(image_path)

    # pad the image
    img = padImage(img)

    # get the contours
    contours = main(img)

    # use contours to create new image...
    im = fillImage(img, contours)
    contour_im = drawContours(img, contours)
    print(f"{contour_im.shape=}")
    tf.imwrite(f"{output_name[:-8]}.tiff", contour_im)
    # save image
    # np.save("control_contours.npy", im)

    idf = img2df(im)
    
    # create 3D polygon...
    final_geom = idf.getPolygon3D()

    gdf = idf.getGDF()

    # save final boundaries...
    gpd.GeoSeries([final_geom]).to_file(f"{output_name}", driver="GeoJSON")

    # save the gdf dataframe
    gdf.to_file(f"DF{output_name}",driver="GeoJSON")
