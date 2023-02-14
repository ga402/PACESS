import numpy as np
import tifffile as tf
import matplotlib.pyplot as plt
import cv2 as cv
import seaborn as sns
from scipy.interpolate import splprep, splev
from utils.func import *
from functools import partial
import functools


@apply3D
def padImage(img):
    """Pads the image with a line of zeros on the 4 edges

    Args:
        img (int): binary image

    Returns:
        img (int): binary image with pad of zeros.
    """
    return np.pad(img, [(1, 1), (1, 1)], mode="constant")






@papply3D
def padImagePar(img):
    """Pads the image with a line of zeros on the 4 edges

    Args:
        img (int): binary image

    Returns:
        img (int): binary image with pad of zeros.
    """
    return np.pad(img, [(1, 1), (1, 1)], mode="constant")







def findContours(binary_img):
    """extracts contours and hierarchy using tree function

    Args:
        binary_img ([np.array]): binary image post padding

    Returns:
        contours: [description]
        hierarchy: [description]
    """
    contours, hierarchy = cv.findContours(
        binary_img, cv.RETR_TREE, cv.CHAIN_APPROX_SIMPLE
    )
    return contours, hierarchy






def getContourAreasPerimeters(contour):
    """get area and perimeter length of all the contours you extract

    Args:
        contour ([type]): [description]

    Returns:
        [type]: [description]
    """
    area = cv.contourArea(contour)
    perimeter = cv.arcLength(contour, True)
    return area, perimeter






def extractAreasPerimeters(contours):
    """apply `getContourAreasPerimeters` to each the contours

    Args:
        contours ([type]): [description]

    Returns:
        areas and perimeters
    """
    area, perimeter = zip(*map(getContourAreasPerimeters, contours))
    return area, perimeter






def binaryArea(contour, a=5000):
    """identify is contour above a threshold or not

    Args:
        contour ([type]): [description]
        a (int, optional): [description]. Defaults to 5000.

    Returns:
        boolean: True/False
    """
    area = cv.contourArea(contour)
    if area > a:
        return True
    else:
        return False






def filterContoursSize(contours, hierarchy, number=10):
    """filter contours by size

    Args:
        contours ([type]): [description]
        hierarchy ([type]): [description]
        number (int, optional): number of contours to return . Defaults to 10.

    Returns:
        filtered contours
        filtered hierarchy
    """
    areas, _ = extractAreasPerimeters(contours)
    # order
    areas = np.array(areas)
    areas.sort()
    # filter
    filtered_contours = [
        c for c in contours if binaryArea(c, areas[-(number - 1) :][0])
    ]
    filtered_hierarchy = [
        h
        for c, h in zip(contours, hierarchy)
        if binaryArea(c, areas[-(number - 1) :][0])
    ]
    #
    return filtered_contours, filtered_hierarchy







def smoothContour(contour, smooth_parameter=10):
    """smooth the contours using interpolation and a slightly reduced length of contours

    Args:
        contour ([numpy.ndarray]): coords
        smooth_parameter (int, optional): [description]. Defaults to 10.

    Returns:
        [type]: smoothed contours
    """
    endlength = contour.shape[0] - contour.shape[0] % smooth_parameter
    x, y = contour.T
    # Convert from numpy arrays to normal arrays
    x = x.tolist()[0]
    y = y.tolist()[0]
    tck, u = splprep([x, y], u=None, s=1.0, per=1, quiet=1)  # interpolation
    u_new = np.linspace(u.min(), u.max(), endlength)
    x_new, y_new = splev(u_new, tck, der=0)
    # back to np.array
    res_array = [[[int(i[0]), int(i[1])]] for i in zip(x_new, y_new)]
    return np.array(res_array)







def applySmoothContour(contours, smooth_parameter=10):
    """application of the smooth contour function to each contour

    Args:
        contours (list): list of contour
        smooth_parameter (int, optional): [description]. Defaults to 10.

    Returns:
        list: smoothed_contours
    """
    smoothed_contours = list(
        map(partial(smoothContour, smooth_parameter=smooth_parameter), contours)
    )
    return smoothed_contours






@applyVec2Arr3D
def fillImage(img, contours):
    """Create a mask of contours

    Args:
        img ([type]): [description]
        contours ([type]): [description]

    Returns:
        [type]: [description]
    """
    stencil = np.zeros(img.shape).astype("uint8")
    color = [255, 255, 255]
    return cv.fillPoly(stencil, contours, color)






@applyVec2Arr3D
def contours2Image(img, contours):
    """Add contours onto image

    Args:
        img ([type]): [description]
        contours ([type]): [description]

    Returns:
        [type]: [description]
    """
    im = cv.cvtColor(img.copy(), cv.COLOR_GRAY2RGB)
    im = cv.drawContours(im, contours, -1, (255,0 ,0), 3)
    return im






@applyVec2Arr3D
def drawContours(img, contours):
    """Draw a linear contours (no fill)

    Args:
        img ([type]): [description]
        contours ([type]): [description]

    Returns:
        [type]: [description]
    """
    stencil = np.zeros(img.shape).astype("uint8")
    color = [255, 255, 255]
    stencil = cv.drawContours(stencil, contours, -1, color, 3)
    return stencil






@sapply3D
def img2Contours(binary_img):
    """extracts contours and hierarchy using tree function

    Args:
        binary_img ([np.array]): binary image post padding

    Returns:
        contours: [description]
        hierarchy: [description]
    """
    contours, _ = cv.findContours(
        binary_img, cv.RETR_TREE, cv.CHAIN_APPROX_SIMPLE
    )
    return contours