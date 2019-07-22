"""Test module for wrapping c++"""
from __future__ import division
cimport cython
from mesh import VesselNetwork
import numpy as np
cimport numpy as np

DTYPE = np.float

ctypedef np.float_t DTYPE_t
ctypedef np.ndarray[DTYPE_t, ndim=1] (*f_mesh_type)(np.ndarray[DTYPE_t, ndim=2], DTYPE,
                                                    DTYPE, DTYPE, int)
# ctypedef VesselNetworkBase MESHTYPE

cdef class _BcFunc:
   cdef np.ndarray[DTYPE_t, ndim=1] (*func)(self, np.ndarray[DTYPE_t, ndim=2],
                                            DTYPE_t, DTYPE_t, DTYPE_t, np.int_t)

# @cython.boundscheck(False) # turn off bounds-checking for entire function
# @cython.wraparound(False)  # turn off negative index wrapping for entire function
def naive_sum(np.ndarray[DTYPE_t, ndim=1] f, np.ndarray[DTYPE_t, ndim=1] g):
    assert f.dtype == DTYPE and g.dtype == DTYPE
    cdef int vmax = f.shape[0]
    cdef int smax = g.shape[0]

    cdef np.ndarray[DTYPE_t, ndim=1] h = np.zeros(vmax, dtype=DTYPE)

    cdef DTYPE_t value
    cdef int x
    for x in range(vmax):
        h[x] = f[x] + g[x]
    return h
