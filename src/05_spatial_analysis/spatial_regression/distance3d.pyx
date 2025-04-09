# cython: boundscheck=False, wraparound=False, cdivision=True
import numpy as np
cimport numpy as np
from libc.math cimport sqrt, pow

def distance3d(np.ndarray[np.double_t, ndim=2] mat):
    cdef int nrows = mat.shape[0]
    cdef int ncols = mat.shape[1]

    if ncols != 3:
        raise ValueError("Incompatible number of dimensions")

    cdef np.ndarray[np.double_t, ndim=2] dmat = np.zeros((nrows, nrows), dtype=np.float64)
    cdef int r1, r2, c
    cdef double total, diff

    for r1 in range(nrows):
        for r2 in range(nrows):
            total = 0.0
            for c in range(3):  # since we know ncols == 3
                diff = mat[r1, c] - mat[r2, c]
                total += diff * diff
            dmat[r1, r2] = sqrt(total)

    return dmat
