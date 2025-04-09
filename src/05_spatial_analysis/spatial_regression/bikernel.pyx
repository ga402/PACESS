# cython: boundscheck=False, wraparound=False, cdivision=True
import numpy as np
cimport numpy as np
from libc.math cimport pow

def bikernel(np.ndarray[np.double_t, ndim=2] dmat, double b):
    cdef int nrows = dmat.shape[0]
    
    if dmat.shape[0] != dmat.shape[1]:
        raise ValueError("dmat should be a square matrix")

    cdef np.ndarray[np.double_t, ndim=2] bikernel_mat = np.zeros((nrows, nrows), dtype=np.float64)
    
    cdef int i, j
    cdef double val, ratio

    for i in range(nrows):
        for j in range(nrows):
            val = dmat[i, j]
            if val < b:
                ratio = val / b
                bikernel_mat[i, j] = pow(1 - pow(ratio, 2), 2)
            else:
                bikernel_mat[i, j] = 0.0

    return bikernel_mat
