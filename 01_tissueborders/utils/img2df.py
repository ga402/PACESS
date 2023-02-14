# type: ignore
import numpy as np
import pandas as pd
import json
import os
import glob
import tifffile as tf
from shapely.geometry import Polygon, MultiPolygon
import geopandas as gpd
import matplotlib.pyplot as plt
import cv2 as cv
from utils.func import *


class baseImage:
    def __init__(self, image):
        try:
            self.image = image
            self.z, self.y, self.x = self.image.shape
            self.z_levels = range(self.z)
            print(
                f"image loaded.\ny={self.image.shape[1]}\nx={self.image.shape[2]}\nz={self.image.shape[0]}"
            )
        except:
            print("unable to load image file")

    def showImage(self, level, save=False, figure_name="figure"):
        if level > self.image.shape[0] or level < 0:
            return print(
                f"level outside the z dimensions of image  (0 to {image.shape[0]})"
            )
        fig, ax = plt.subplots(figsize=(10, 10))
        ax.imshow(self.image[level, :, :], cmap="gray")
        plt.axis("off")
        if save == True:
            plt.savefig(f"{figure_name}.png")
        plt.show()

    def getImage(self):
        return self.image

    def printShape(self):
        print(
            f"\ny={self.image.shape[1]}\nx={self.image.shape[2]}\nz={self.image.shape[0]}"
        )

    def getDim(self):
        return self.z, self.y, self.x


class img2df(baseImage):
    def __init__(self, image):
        super().__init__(image)
        self.final_geom = None
        self.gdf = None

    @staticmethod
    def _extractCont(im):
        cont, heir = cv.findContours(im, cv.RETR_TREE, cv.CHAIN_APPROX_SIMPLE)
        return cont, heir

    @staticmethod
    def fPolygon(cont):
        cont = cont.squeeze()
        return Polygon(cont)

    def _convertToPolygons(self, contours):
        polygons = map(self.fPolygon, contours)  # converting to Polygons
        return MultiPolygon(polygons)

    def _contours2Polygons(self, im):
        c, h = self._extractCont(im)
        return self._convertToPolygons(c)

    def contours2Polygons(self, image):
        assert len(image.shape) == 3, "wrong dimensions"
        MP = [self._contours2Polygons(image[i, ...]) for i in range(image.shape[0])]
        return MP

    def generateContoursDataFrame(self):
        MP_list = self.contours2Polygons(self.image)
        self.gdf = gpd.GeoDataFrame({"z_level": self.z_levels, "geometry": MP_list})
        self.gdf.set_geometry("geometry")
        print(self.gdf.head(5))

    def getGDF(self):
        if self.gdf is None:
            self.generateContoursDataFrame()
        # self.gdf.set_geometry("geometry")
        return self.gdf.set_geometry("geometry")

    def generateCompleteMultiPolygon(self):
        """Create a single multiple polygon (at each z level)"""
        final_geom = []
        for g, z in zip(self.gdf["geometry"], self.gdf["z_level"]):
            for polygon in g:
                p = Polygon([t + (float(z),) for t in list(polygon.exterior.coords)])
                final_geom.append(p)
            # final_geom.append(MultiPolygon(empty_p))
        final_geom = MultiPolygon(final_geom)
        res = np.where(final_geom.has_z, "has z axis", " has no z axis")
        res = str(res)
        print(f"{final_geom.geom_type} generated {res}")
        self.final_geom = final_geom

    def getPolygon3D(self):

        # if any(self.gdf == None):
        try:
            self.generateContoursDataFrame()
        except:
            print("generateContoursDataFrame failed")

        # if self.final_geom == None:
        try:
            self.generateCompleteMultiPolygon()
        except:
            print("generateCompleteMultiPolygon failed")

        return self.final_geom
