import pandas as pd
import geopandas as gpd
import numpy as np
import glob as glob
from itertools import starmap
from functools import partial
import argparse
import sys
from pathos.multiprocessing import ProcessPool


def get_program_parameters():
    epilogue = """
    create a dataframe of pixel sizes and z levels and image names
    """
    parser = argparse.ArgumentParser(
        description="create the dataframe",
        epilog=epilogue,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-f",
        "--files",
        nargs="+",
        default=["file1.csv", "file2.csv"],
        help="input files",
    )
    parser.add_argument(
        "-c",
        "--cells",
        nargs="+",
        default=["AML", "C", "MGK", "O", "TC"],
        help="input cells",
    )
    parser.add_argument(
        "--xyz_correction",
        nargs="+",
        type=float,
        default=[0.41, 0.41, 5],
        help="input files",
    )
    parser.add_argument("-e", "--edge_length", type=int, help="radius of point freq")
    parser.add_argument(
        "-g", "--gdf", help="geopandas dataframe", default="file.geoJson"
    )
    parser.add_argument("-o", "--output", help="output csv name, eg. output.csv")
    parser.add_argument(
        "-d", "--data_out", help="name of the output filtered data, eg out.csv"
    )
    parser.add_argument("--save_data_out", type=str)
    args = parser.parse_args()
    return (
        args.files,
        args.cells,
        args.xyz_correction,
        args.edge_length,
        args.gdf,
        args.output,
        args.data_out,
        args.save_data_out,
    )


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ("True", "TRUE", "yes", "true", "t", "y", "1"):
        return True
    elif v.lower() in ("False", "FALSE", "no", "false", "f", "n", "0"):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")


def readGDF(gdf_path):
    if not gdf_path.endswith(".geojson"):
        raise Exception("geopandas dataframe required with .geojson")
    try:
        gdf = gpd.read_file(gdf_path)
        assert isinstance(gdf, gpd.GeoDataFrame), "no a geopandas dataframe"
    except:
        print("Unable to load geopandas dataframe")
        sys.exit()
    return gdf


def readFile(data_path, label):
    try:
        df = pd.read_csv(data_path)
        df["cell"] = label
    except:
        print("Unable to load file")
        sys.exit()
    return df


def generatePointsFromXY(xints, yints, index, geom):
    gdf_points = gpd.GeoSeries(gpd.points_from_xy(xints, yints))
    points = gpd.GeoDataFrame({"geometry": gdf_points, "Index": index})
    poly = gpd.GeoDataFrame({"geometry": geom})
    within_points = gpd.sjoin(points, poly, op="within")
    #within_points = within_points.reset_index().drop("index", axis=1)
    return within_points



def removeOutsideLattice(gdf_row, xyz, edge_length):
    assert xyz.shape[1] == 3, "'xyz needs to be an nx3 array"
    assert isinstance(gdf, gpd.GeoDataFrame), "no a geopandas dataframe"
    assert "z_level" in gdf.columns, "geopandas needs column called 'z_level'"
    z = gdf_row["z_level"]
    xy = xyz[np.where(xyz[:, 2] == z)]
    geom_geometry = gdf_row["geometry"].convex_hull.buffer(edge_length)
    within_points_xy = generatePointsFromXY(xy[:, 0], xy[:, 1], z, [geom_geometry])
    return within_points_xy


def generateLattice(xmin, xmax, ymin, ymax, zmax, edge_length, z_correction=5):
    # assert edge_length < zmax*z_correction, "'r' should be smaller than the max z value"
    assert edge_length < xmax - xmin, "'r' should be smaller than the max x value"
    assert edge_length < ymax - ymin, "'r' should be smaller than the max y value"
    x_ = np.arange(xmin, xmax, edge_length)
    y_ = np.arange(ymin, ymax, edge_length)
    z_ = np.arange(0, zmax, 1)  ## going to take every z layer
    x, y, z = np.meshgrid(x_, y_, z_, indexing="ij")
    xyz = list(map(lambda x: x[:, np.newaxis].ravel(), [x, y, z]))
    return np.c_[xyz].T


def removeOutsideDataFrame(gdf_row, df):
    z = gdf_row["z_level"]
    dfZ = df.query(f"z_approx == {z}")
    geom = gdf_row["geometry"]
    within_points = generatePointsFromXY(dfZ.x, dfZ.y, dfZ.index, geom)
    return dfZ.loc[within_points.Index]



def getBoundsLattice(gdf):
    zmax = gdf.__len__()
    bounds_arr = np.array([[gdf.geometry[i].bounds] for i in range(zmax)]).squeeze()
    xmin, ymin = tuple(np.min(bounds_arr[:, [0, 1]], axis=0))
    xmax, ymax = tuple(np.max(bounds_arr[:, [2, 3]], axis=0))
    return xmin, xmax, ymin, ymax, zmax


def getXYZ(lattice):
    x = lattice["geometry"].x
    y = lattice["geometry"].y
    return np.c_[x, y, lattice["Index"]]


def getBoxLimits(lattice_corrected, radius, xyz_correction): 
    radius_correction = radius*xyz_correction[0]
    xl_min, xl_max = lattice_corrected[:, 0]-radius_correction, lattice_corrected[:, 0]+radius_correction
    yl_min, yl_max = lattice_corrected[:, 1]-radius_correction, lattice_corrected[:, 1]+radius_correction
    zl_min, zl_max = lattice_corrected[:, 2], lattice_corrected[:, 2]+xyz_correction[2]
    return np.c_[xl_min, xl_max, yl_min, yl_max, zl_min, zl_max]


def getNumber(boxlimits, mat):
    return mat.query(f'X >= {boxlimits[0]} & X < {boxlimits[1]} & Y >= {boxlimits[2]} & Y < {boxlimits[3]} & Z >= {boxlimits[4]} & Z < {boxlimits[5]}').shape[0]




def main_function(data_path, gdf, labels, xyz_correction, edge_length, prune=True):
    """_summary_

    Args:
        data_path (_type_): _description_
        labels (_type_): _description_
        xyz (_type_): _description_
        N (_type_): _description_
        radius (_type_): _description_

    Returns:
        _type_: _description_
    """

    args_input = list(zip(data_path, labels))
    dat = pd.concat(list(starmap(readFile, args_input)))
    radius = edge_length / 2
    #
    try:
        dat["z_approx"] = round(dat.zmean)
    except:
        dat["z_approx"] = np.ceil(np.floor(dat.zmean))
    #
    # check asserts
    assert gdf.z_level.max() == dat.zmax.max(), "max z levels don't match"
    assert gdf.z_level.min() == dat.zmin.min(), "min z levels don't match"
    #
    try:
        df = gdf.apply(removeOutsideDataFrame, df=dat, axis=1)
    except:
        raise Exception("unable to filter using polygons")
    #
    # create a simplified mat of the data.
    #
    df = pd.concat(list(df))
    mat = df.loc[:, ["x", "y", "zmean"]] * xyz_correction
    mat = mat.rename(columns={"x": "X", "y": "Y", "zmean": "Z"})
    mat["cell"] = df.cell
    #
    # generate full lattice
    z_c = xyz_correction[2]  # z_correction
    xyz = generateLattice(*getBoundsLattice(gdf), edge_length, z_correction=z_c)
    # xyz[:, 2].max()
    # xyz = generateLattice(2000, 2500, 2500, 3000, 49, edge_length, z_correction=z_c)
    # xyz.shape
    #
    if prune:
        try:
            lattice = gdf.apply(removeOutsideLattice, xyz=xyz, edge_length=edge_length, axis=1)
            lattice = pd.concat(list(lattice)).reset_index().drop("index", axis=1)
            lattice = getXYZ(lattice)  # to array
        except:
            raise Exception("unable to create prune lattice..")
    else:
        lattice = xyz
    #
    print(f"{lattice.shape=}")
    # correct for dimensions
    assert lattice.shape[1] == len(xyz_correction), "lattice shape must be 3"
    lattice_corrected = lattice * xyz_correction  # correction for distance
    #
    boxlimits = getBoxLimits(lattice_corrected, radius, xyz_correction)
    # generative final dataset
    lattice = pd.DataFrame(lattice, columns=["x", "y", "z"])  # final dataframe
    # main loop
    for c in labels:
        subset_df = mat.query(f"cell == '{c}'") ## corrected size
        pool = ProcessPool(nodes=4)
        frequency_values = pool.map(partial(getNumber, mat=subset_df), boxlimits)# adjust to the actual area
        lattice[f"Freq_{c}"] = frequency_values
    #
    return lattice, df


if __name__ == "__main__":
    print("starting...")
    (
        data_path,
        labels,
        xyz_correction,
        edge_length,
        gdf_path,
        output_name,
        df_out,
        save_df,
    ) = get_program_parameters()

    save = str2bool(save_df)
    print(f"Save DataFrame= {save}")
    assert len(data_path) == len(labels), "number of files must equal number of cells"
    gdf = readGDF(gdf_path)
    lattice, df = main_function(data_path, gdf, labels, xyz_correction, edge_length, prune=True)
    # cut down the size of the dataframe by removing all the areas without any cells...
    lattice['cell_n'] = lattice.drop(['x', 'y', 'z'], axis=1).apply(sum, axis=1)
    lattice = lattice.query('cell_n >0').reset_index().drop(['index', 'cell_n'], axis=1)
    # 
    try:
        assert lattice.Freq_AML.sum() == df.query("cell=='AML'").shape[0], "Freq_AML does not equal AML numbers"
        assert lattice.Freq_MGK.sum() == df.query("cell=='MGK'").shape[0],  "Freq_MGK does not equal MGK numbers"
        assert lattice.Freq_TC.sum() == df.query("cell=='TC'").shape[0],  "Freq_TC does not equal TC numbers"
        print("Frequencires of AML, MGK and TC equal counts")
    except:
        print("Frequencies fail to match counts")
        xt = df.query("cell=='AML'").shape[0]
        print(f" >  estimated AML: {lattice.Freq_AML.sum()} - True AML :{xt}")
        xt = df.query("cell=='MGK'").shape[0]
        print(f" >  estimated MGK: {lattice.Freq_MGK.sum()} - True MGK :{xt}")
        xt = df.query("cell=='TC'").shape[0]
        print(f" >  estimated TC: {lattice.Freq_TC.sum()} - True TC :{xt}")
    if save:
        df.to_csv(f"{df_out}", index=False)
    print(f"Final lattice size = {lattice.shape}")
    lattice.to_csv(f"{output_name}", index=False)
    print(f"complete lattice with edge_length = {edge_length}")
    sys.exit()



