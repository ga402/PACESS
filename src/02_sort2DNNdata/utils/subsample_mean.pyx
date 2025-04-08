# distutils: language = c
import numpy as np
cimport numpy as cnp
cimport cython
from cython.parallel import prange

cnp.import_array()  # Initialize NumPy support

@cython.boundscheck(False)  # Disable bounds checking for speed
@cython.wraparound(False)   # Disable negative indexing overhead
@cython.cdivision(True)     # Enable C-style division for speed
def subsample_mean(cnp.ndarray[signed char, ndim=2] arr, int x_start, int y_start, int width, int height):
    """
    Subsamples a 2D NumPy array and returns the mean value of the selected region.

    Args:
        arr (np.ndarray): Input 2D NumPy array.
        x_start (int): Starting x-coordinate (row index).
        y_start (int): Starting y-coordinate (column index).
        width (int): Width of the sub-region.
        height (int): Height of the sub-region.

    Returns:
        float: Mean value of the selected region.
    """
    cdef int x, y
    cdef double sum_val = 0.0
    cdef int count = 0

    for x in range(x_start, x_start + width):
        for y in range(y_start, y_start + height):
            sum_val += arr[x, y]
            count += 1
    
    return sum_val / count if count > 0 else 0.0 



@cython.boundscheck(False)  # Disable bounds checking for performance
@cython.wraparound(False)   # Disable negative indexing
@cython.cdivision(True)     # Enable C-style division
def rgb_to_grayscale(cnp.ndarray[signed char, ndim=3] rgb_image):
    """
    Converts an RGB image (H, W, 3) to a Grayscale image (H, W) using the luminance formula.

    Args:
        rgb_image (np.ndarray): 3D NumPy array (H, W, 3) with RGB channels.

    Returns:
        np.ndarray: 2D NumPy array (H, W) representing the grayscale image.
    """
    cdef int height = rgb_image.shape[0]
    cdef int width = rgb_image.shape[1]
    
    # Create grayscale output array
    cdef cnp.ndarray[signed char, ndim=2] gray_image = np.empty((height, width), dtype=np.int8)

    cdef int i, j
    cdef double r, g, b, gray_value

    # Convert RGB to grayscale using the luminance formula
    for i in range(height):
        for j in range(width):
            r = rgb_image[i, j, 0]
            g = rgb_image[i, j, 1]
            b = rgb_image[i, j, 2]

            # Compute grayscale value using standard weights
            gray_value = 0.2989 * r + 0.5870 * g + 0.1140 * b
            
            # Assign to output (clamp values between 0-255)
            gray_image[i, j] = <signed char> gray_value

    return gray_image




def vectorized_calc_mean_from_image_dict(dict image_dict, cnp.ndarray[object, ndim=1] image_names, cnp.ndarray[cnp.uint64_t, ndim=1] x_pos, cnp.ndarray[cnp.uint64_t, ndim=1] y_pos, cnp.ndarray[cnp.uint64_t, ndim=1] width, cnp.ndarray[cnp.uint64_t, ndim=1] height, int color):
    """
    Highly optimized, parallelized function applied to a NumPy array.

    Args:
        arr (np.ndarray): 1D NumPy array.

    Returns:
        np.ndarray: Transformed NumPy array.
    """
    cdef Py_ssize_t i, size = x_pos.shape[0]
    cdef cnp.ndarray[cnp.float64_t, ndim=1] result = np.empty(size, dtype=np.float64)
    cdef cnp.ndarray[signed char, ndim=3] img_array
    cdef cnp.ndarray[signed char, ndim=2] img_col
    # 
    
    for i in range(size): 
        img_array = np.ascontiguousarray(image_dict[image_names[i]], dtype=np.int8) 

        if color > 3:
            img_col = rgb_to_grayscale(img_array)
        else:
            img_col = img_array[:, :, color]

        result[i] = subsample_mean(img_col, x_pos[i], y_pos[i], width[i], height[i])


    return result