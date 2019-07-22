cimport cython
import numpy as np
cimport numpy as np
cimport cymesh
cimport c_cynum
from c_cynum cimport Vessel
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.map cimport map
from cython.operator cimport dereference as deref


cdef class cyVessel:
    # cdef cymesh.Vessel_t *thisptr
    def __cinit__(self, string name, double L, double R_prox, double R_dist,
                  double Wall_th, dict Windk = dict(), int Id = 0, double rc=0):
        if type(self) is cyVessel:
            self.thisptr = new cymesh.Vessel_t(name, L, R_prox, R_dist, Wall_th, Windk, Id)

    cpdef name(self):
        return self.thisptr.getName()

    def __dealloc__(self):
        if type(self) is cyVessel:
            del self.thisptr
