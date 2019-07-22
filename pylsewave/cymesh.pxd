cimport cython
cimport c_cynum
from c_cynum cimport Vessel
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.map cimport map

# cdef class cyVessel(object)

ctypedef c_cynum.Vessel Vessel_t

cdef class cyVessel(object):
    cdef Vessel_t *thisptr

    cpdef name(self)
    # cpdef __cinit__(self, string name, double L, double R_prox, double R_dist,
    #               double Wall_th, dict Windk, int Id, double rc)
    # cpdef __dealloc__(self)