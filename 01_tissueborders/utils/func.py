import numpy as np
import sys
import argparse
from datetime import datetime
import logging
import os
import cv2 as cv
import functools
from PIL import Image
import matplotlib.pyplot as plt

## for parallel proccessing
from joblib import delayed, Parallel
import multiprocessing


def createlog():
    fn = __file__.replace(".py", "")
    logfile = datetime.now().strftime("logfile%Hhr%Mmin%d%m%Y.log")
    logfile = fn + logfile
    return logfile


def dfuncapply(func, img):
    return np.array([func(img[i, ...]) for i in np.arange(0, img.shape[0])])


def sfuncapply(func, img):
    return [func(img[i, ...]) for i in np.arange(0, img.shape[0])]


def dfuncapply2Arr(func, img, img2):
    return np.array(
        [func(img[i, ...], img2[i, ...]) for i in np.arange(0, img.shape[0])]
    )


def dfuncapplyVec2Arr(func, img, vec):
    assert img.shape[0] == len(vec), "array and vector dimensions do not match"
    return np.array([func(img[i, ...], vec[i]) for i in np.arange(0, img.shape[0])])


def apply3D(func):
    """function to allow you to apply function in 3D"""

    @functools.wraps(func)
    def wrapper_apply3D(*args, **kwargs):
        return dfuncapply(func, *args)

    return wrapper_apply3D


def sapply3D(func):
    """function to allow you to apply function in 3D"""

    @functools.wraps(func)
    def wrapper_sapply3D(*args, **kwargs):
        return sfuncapply(func, *args)

    return wrapper_sapply3D


def apply2Arr3D(func):
    """function to allow you to apply function in 3D"""

    @functools.wraps(func)
    def wrapper_apply3D2Arr(*args, **kwargs):
        return dfuncapply2Arr(func, *args)

    return wrapper_apply3D2Arr


def applyVec2Arr3D(func):
    """function to allow you to apply function in 3D"""

    @functools.wraps(func)
    def wrapper_apply3DVec2Arr(*args, **kwargs):
        return dfuncapplyVec2Arr(func, *args)

    return wrapper_apply3DVec2Arr


@apply3D
def createMask(img):
    """
    This function converts the 4-dimension image to the 3D image.
    """
    blur = cv.GaussianBlur(img.copy(), (5, 5), 0)
    _, bin = cv.threshold(blur, 0, 255, cv.THRESH_BINARY + cv.THRESH_OTSU)
    return bin


def getCPUs():
    """get the number of CPUs available

    Returns:
        [int]: number of CPUs on machine
    """
    try:
        cpus = multiprocessing.cpu_count()
        cpus = int(cpus / 2)
    except NotImplementedError:
        cpus = 2  # arbitrary default
    finally:
        return cpus


def parallelfuncapply(func, img):
    """apply function in parallel

    Args:
        func ([type]): [description]
        img ([type]): [description]

    Returns:
        [type]: [description]
    """
    cpu_count = getCPUs()
    values = Parallel(n_jobs=cpu_count)(
        delayed(func)(img[i, ...]) for i in np.arange(0, img.shape[0])
    )
    return np.array(values)


def sparallelfuncapply(func, img):
    """apply function in parallel

    Args:
        func ([type]): [description]
        img ([type]): [description]

    Returns:
        [type]: [description]
    """
    cpu_count = getCPUs()
    values = Parallel(n_jobs=cpu_count)(
        delayed(func)(img[i, ...]) for i in np.arange(0, img.shape[0])
    )
    return values


def papply3D(func):
    """decorator of parallelfuncapply

    Args:
        func ([type]): [description]

    Returns:
        [type]: [description]
    """

    @functools.wraps(func)
    def wrapper_papply3D(*args, **kwargs):
        return parallelfuncapply(func, *args)

    return wrapper_papply3D


def spapply3D(func):
    """decorator of parallelfuncapply

    Args:
        func ([type]): [description]

    Returns:
        [type]: [description]
    """

    @functools.wraps(func)
    def wrapper_spapply3D(*args, **kwargs):
        return sparallelfuncapply(func, *args)

    return wrapper_spapply3D
