
import numpy as np
import pandas as pd
from matplotlib.text import Annotation
from mpl_toolkits.mplot3d.axes3d import Axes3D
from mpl_toolkits.mplot3d.proj3d import proj_transform
import math
import scipy
from scipy.signal import fftconvolve




# functions
def checkxyratio(dt):
    xx, yy, _ = dt.loc[:, ['x','y','z']].apply(lambda x: abs(x.max()-  x.min()))
    if yy == 0.0:
        return 0.0001
    elif xx == 0.0:
        return 0.0001
    else:
        return xx/yy

def sort_order_clusterSize(df):
    df['cluster_order'] = df.cluster.map(df.groupby('cluster').count().x.to_dict())
    df = df.sort_values('cluster_order', ascending=False)
    return df

class Annotation3D(Annotation):
    def __init__(self, text, xyz, *args, **kwargs):
        super().__init__(text, xy=(0, 0), *args, **kwargs)
        self._xyz = xyz
    def draw(self, renderer):
        x2, y2, z2 = proj_transform(*self._xyz, self.axes.M)
        self.xy = (x2, y2)
        super().draw(renderer)


def _annotate3D(ax, text, xyz, *args, **kwargs):
    '''Add anotation `text` to an `Axes3d` instance.'''
    annotation = Annotation3D(text, xyz, *args, **kwargs)
    ax.add_artist(annotation)

setattr(Axes3D, 'annotate3D', _annotate3D)


def roundup(x):
    return int(math.ceil(x / 100.0)) * 100


def filterValues(x, y, v, gdf, z_level):
    import geopandas as gpd
    assert isinstance(gdf, gpd.GeoDataFrame),f"{gdf} needs to be geopandas df"
    temp = pd.DataFrame({"x":x, "y":y, "v": v})
    gdf_points = gpd.GeoSeries(gpd.points_from_xy(temp.x.values, temp.y.values))
    points = gpd.GeoDataFrame({"geometry": gdf_points}).set_crs(epsg=4326, inplace=True)
    geom = gdf.query(f"z_level=={z_level}").geometry
    poly = gpd.GeoDataFrame({"geometry": geom}).set_crs(epsg=4326, inplace=True)
    within_points = gpd.sjoin(points, poly, op="within")
    #temp.loc[temp.index.isin(within_points.index),]
    return temp.loc[temp.index.isin(within_points.index),].values



def calculate_hull(
        X, 
        scale=1.1, 
        padding="scale", 
        n_interpolate=100, 
        interpolation="quadratic_periodic", 
        return_hull_points=False):
    """
    Calculates a "smooth" hull around given points in `X`.
    The different settings have different drawbacks but the given defaults work reasonably well.
    Parameters
    ----------
    X : np.ndarray
        2d-array with 2 columns and `n` rows
    scale : float, optional
        padding strength, by default 1.1
    padding : str, optional
        padding mode, by default "scale"
    n_interpolate : int, optional
        number of interpolation points, by default 100
    interpolation : str or callable(ix,iy,x), optional
        interpolation mode, by default "quadratic_periodic"

    Inspired by: https://stackoverflow.com/a/17557853/991496
    """
    if padding == "scale":   
        # scaling based padding
        scaler = sklearn.pipeline.make_pipeline(
            sklearn.preprocessing.StandardScaler(with_std=False),
            sklearn.preprocessing.MinMaxScaler(feature_range=(-1,1)))
        points_scaled = scaler.fit_transform(X) * scale
        hull_scaled = scipy.spatial.ConvexHull(points_scaled, incremental=True)
        hull_points_scaled = points_scaled[hull_scaled.vertices]
        hull_points = scaler.inverse_transform(hull_points_scaled)
        hull_points = np.concatenate([hull_points, hull_points[:1]])
    #
    elif padding == "extend" or isinstance(padding, (float, int)):
        # extension based padding
        # TODO: remove?
        if padding == "extend":
            add = (scale - 1) * np.max([
                X[:,0].max() - X[:,0].min(), 
                X[:,1].max() - X[:,1].min()])
        else:
            add = padding
        points_added = np.concatenate([
            X + [0,add], 
            X - [0,add], 
            X + [add, 0], 
            X - [add, 0]])
        hull = scipy.spatial.ConvexHull(points_added)
        hull_points = points_added[hull.vertices]
        hull_points = np.concatenate([hull_points, hull_points[:1]])
    else:
        raise ValueError(f"Unknown padding mode: {padding}")
    # number of interpolated points
    nt = np.linspace(0, 1, n_interpolate)
    x, y = hull_points[:,0], hull_points[:,1]
    # ensures the same spacing of points between all hull points
    t = np.zeros(x.shape)
    t[1:] = np.sqrt((x[1:] - x[:-1])**2 + (y[1:] - y[:-1])**2)
    t = np.cumsum(t)
    t /= t[-1]
    # interpolation types
    if interpolation is None or interpolation == "linear":
        x2 = scipy.interpolate.interp1d(t, x, kind="linear")(nt)
        y2 = scipy.interpolate.interp1d(t, y, kind="linear")(nt)
    elif interpolation == "quadratic":
        x2 = scipy.interpolate.interp1d(t, x, kind="quadratic")(nt)
        y2 = scipy.interpolate.interp1d(t, y, kind="quadratic")(nt)
    #
    elif interpolation == "quadratic_periodic":
        x2 = scipy.interpolate.splev(nt, scipy.interpolate.splrep(t, x, per=True, k=4))
        y2 = scipy.interpolate.splev(nt, scipy.interpolate.splrep(t, y, per=True, k=4))
    #
    elif interpolation == "cubic":
        x2 = scipy.interpolate.CubicSpline(t, x, bc_type="periodic")(nt)
        y2 = scipy.interpolate.CubicSpline(t, y, bc_type="periodic")(nt)
    else:
        x2 = interpolation(t, x, nt)
        y2 = interpolation(t, y, nt)
    #
    X_hull = np.concatenate([x2.reshape(-1,1), y2.reshape(-1,1)], axis=1)
    if return_hull_points:
        return X_hull, hull_points
    else:
        return X_hull




def filterValues(x, y, v, gdf, z_level):
    import geopandas as gpd
    assert isinstance(gdf, gpd.GeoDataFrame),f"{gdf} needs to be geopandas df"
    temp = pd.DataFrame({"x":x, "y":y, "v": v})
    gdf_points = gpd.GeoSeries(gpd.points_from_xy(temp.x.values, temp.y.values))
    points = gpd.GeoDataFrame({"geometry": gdf_points}).set_crs(epsg=4326, inplace=True)
    geom = gdf.query(f"z_level=={z_level}").geometry
    poly = gpd.GeoDataFrame({"geometry": geom}).set_crs(epsg=4326, inplace=True)
    within_points = gpd.sjoin(points, poly, op="within")
    #temp.loc[temp.index.isin(within_points.index),]
    return temp.loc[temp.index.isin(within_points.index),].values


def getBlurArray(df, z, group, size=1):
    array3  = convertDfToArr(df.query(f'z == {z}'), group).squeeze()
    array4 = np.where(np.isnan(array3), 0, array3)
    array4 = gaussian_blur(array4, size)
    array4 = np.where(np.isnan(array3), np.nan, array4)
    return array4





def gaussian_blur(in_array, size):
    # expand in_array to fit edge of kernel
    padded_array = np.pad(in_array, size, 'symmetric')
    # build kernel
    x, y = np.mgrid[-size:size + 1, -size:size + 1]
    g = np.exp(-(x**2 / float(size) + y**2 / float(size)))
    g = (g / g.sum()).astype(in_array.dtype)
    # do the Gaussian blur
    return fftconvolve(padded_array, g, mode='valid')


import xarray as xr


def convertDfToArr(df, VALUE):
    dt = pd.pivot_table(df, values=f'{VALUE}', index=['x', 'y', 'z'])
    xrTensor = xr.DataArray(dt).unstack("dim_0")
    array = xrTensor.values[0].T
    return array