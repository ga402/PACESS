# distutils: language = c
from PIL import Image
from pathlib import Path
cimport cython
import numpy as np
cimport numpy as cnp 
from libc.stdlib cimport malloc, free
from libc.string cimport memcpy
cnp.import_array() 


@cython.boundscheck(False)  # Disable bounds checking for speed
@cython.wraparound(False)   # Disable negative index wraparound for speed
def load_image_as_numpy(str image_path):
    """
    Loads an image using Pillow and converts it to a NumPy array.

    Args:
        image_path (str): Path to the image file.

    Returns:
        numpy.ndarray: NumPy array representation of the image.
    """
    img = Image.open(image_path).convert("RGB")  # Load image
    cdef cnp.ndarray[cnp.uint8_t, ndim=3] img_array = np.array(img, dtype=np.uint8)  # Convert to NumPy array
    return img_array



def load_images_to_dict(list image_paths):
    """
    Loads images using Pillow, converts to NumPy arrays, and stores them in a dictionary.

    Args:
        image_paths (list): List of image file paths.

    Returns:
        dict: Dictionary with {filename: numpy array}
    """
    cdef dict image_dict = {}
    cdef str img_path, img_name
    #cdef Image.Image img
    cdef cnp.ndarray[cnp.uint8_t, ndim=3] img_array

    for img_path in image_paths:
        # Extract filename
        img_name = str(Path(img_path).stem)
        
        img_array = load_image_as_numpy(img_path)
        
        # Store in dictionary
        image_dict[img_name] = img_array

    return image_dict