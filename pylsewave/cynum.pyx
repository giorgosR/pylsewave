"""
This is a Cython file with classes for pylsewave toolkit. Cython translates everything to C++.
This module contains solvers that support OPENMP (via Cython OpenMP) for parallel CPU calcs.
"""
cimport cython
from cython cimport view
from cython.parallel import prange
from cython.parallel import parallel
# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np
cimport c_cynum
from c_cynum cimport (Vessel, VesselScaled,
                      VesselNetwork, VesselNetworkSc, write2file,
                      tdma, std_dev, grad)
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.map cimport map
from libc.math cimport sqrt, pow, M_PI, fabs, exp
from libc.stdio cimport printf
from libc.stdlib cimport calloc, malloc, free
from cython.operator cimport dereference as deref
import os
from cymem.cymem cimport Pool

__author__ = "Georgios E. Ragkousis"

ctypedef np.float_t DTYPE_t
ctypedef np.int64_t ITYPE_t

cdef ITYPE_t STATUS_OK = 0
cdef ITYPE_t STATUS_ERROR = -1

PRINT_STATUS = False
WRITE_STATUS = True


cdef class cyVessel(object):
    cdef Vessel *thisptr

    def __cinit__(self, string name, DTYPE_t L, DTYPE_t R_prox, DTYPE_t R_dist,
                  DTYPE_t Wall_th, dict Windk = dict(), int Id = 0, rc=0):
        if type(self) is cyVessel:
            self.thisptr = new Vessel(name, L, R_prox, R_dist, Wall_th, Windk, Id)

    def __dealloc__(self):
        if type(self) is cyVessel:
            del self.thisptr

    @property
    def name(self):
        return self.thisptr.getName()

    @property
    def L(self):
        return self.thisptr.getL()

    @property
    def r_prox(self):
        return self.thisptr.getRadius_prox()

    @property
    def r_dist(self):
        return self.thisptr.getRadius_dist()

    @property
    def id(self):
        return self.thisptr.getId()

    @property
    def w_th(self):
        return self.thisptr.getWall_th()

    property RLC:
        def __get__(self):
            return self.thisptr.getRLC()
        def __set__(self, dict dinput):
            self.thisptr.setRLC(dinput)

    property dx:
        def __get__(self):
            return self.thisptr.getdx()
        def __set__(self, DTYPE_t dinput):
            self.thisptr.setdx(dinput)

    @property
    def x(self):
        return np.asarray(self.thisptr.get_x())

    @property
    def r0(self):
        return np.asarray(self.thisptr.getR0())

    @property
    def f_r0(self):
        return np.asarray(self.thisptr.get_f_R0())

    @property
    def df_dr0(self):
        return np.asarray(self.thisptr.get_df_dR0())

    @property
    def df_dx(self):
        return np.asarray(self.thisptr.get_df_dx())

    @property
    def f_r0_ph(self):
        return np.asarray(self.thisptr.get_f_R0_ph())

    @property
    def df_dr0_ph(self):
        return np.asarray(self.thisptr.get_df_dR0_ph())

    @property
    def df_dx_ph(self):
        return np.asarray(self.thisptr.get_df_dx_ph())

    @property
    def f_r0_mh(self):
        return np.asarray(self.thisptr.get_f_R0_mh())

    @property
    def df_dr0_mh(self):
        return np.asarray(self.thisptr.get_df_dR0_mh())

    @property
    def df_dx_mh(self):
        return np.asarray(self.thisptr.get_df_dx_mh())

    @property
    def k(self):
        return np.asarray(self.thisptr.get_k_vector())

    def set_k_vector(self, input_v):
        self.thisptr.set_k_vector(input_v)

    def interpolate_R0(self, value):
        return np.asarray(self.thisptr.interpolate_R0(value))


cdef class cyVesselScaled(cyVessel):
    cdef VesselScaled *thisptrDerived
    def __cinit__(self, string name, DTYPE_t L, DTYPE_t R_prox, DTYPE_t R_dist,
                  DTYPE_t Wall_th, dict Windk=dict(), int Id=0, DTYPE_t rc=1.0):
        if type(self) is cyVesselScaled:
            self.thisptrDerived = self.thisptr = new VesselScaled(name, L, R_prox, R_dist, Wall_th, Windk, Id, rc)

    def set_k_vector(self, input_v, rho=0, qc=0, rc=0):
        self.thisptrDerived.set_k_vector(input_v, rho, qc, rc)

    def __dealloc__(self):
        if type(self) is cyVesselScaled:
            del self.thisptrDerived

cdef class cyVesselNetwork(object):
    cdef VesselNetwork *thisptr
    cdef object _vessels
    def __cinit__(self, cyVessel vec = None, list vessels=None, DTYPE_t p0=-1, DTYPE_t rho=-1,
                  DTYPE_t Re=-1, DTYPE_t dx=-1, DTYPE_t Nx=-1, DTYPE_t qc=1., DTYPE_t rc=1.):
        cdef vector[Vessel] v = vector[Vessel]()
        cdef cyVessel a
        if type(self) is cyVesselNetwork:
            if vec is not None and vessels is None:
                self.thisptr = new VesselNetwork(deref(vec.thisptr), p0, rho, Re, dx, Nx)
            if vessels is not None and vec is None:
                self._vessels = vessels
                for a in vessels:
                    v.push_back(deref(a.thisptr))
                self.thisptr = new VesselNetwork(v, p0, rho, Re, dx, Nx)

    @property
    def vessels(self):
        return self._vessels

    @property
    def dx(self):
        return self.thisptr.get_dx()

    @dx.setter
    def dx(self, indx):
        self.thisptr.set_dx(indx)

    @property
    def p0(self):
        return self.thisptr.get_p0()

    @p0.setter
    def p0(self, inp0):
        self.thisptr.set_p0(inp0)

    @property
    def Re(self):
        return self.thisptr.get_Re()

    @Re.setter
    def Re(self, inRe):
        self.thisptr.set_Re(inRe)

    @property
    def rho(self):
        return self.thisptr.get_rho()

    @rho.setter
    def rho(self, inrho):
        self.thisptr.set_rho(inrho)

    @property
    def delta(self):
        return self.thisptr.get_delta()

    @delta.setter
    def delta(self, indelta):
        self.thisptr.set_delta(indelta)

    def __dealloc__(self):
        if type(self) is cyVesselNetwork:
            del self.thisptr

cdef class cyVesselNetworkSc(cyVesselNetwork):
    cdef VesselNetworkSc *thisptrDerived
    def __cinit__(self, cyVessel vec = None, list vessels=None, DTYPE_t p0=-1, DTYPE_t rho=-1,
                 DTYPE_t Re=-1, DTYPE_t dx=-1, DTYPE_t Nx=-1, DTYPE_t qc=1., DTYPE_t rc=1.):

        cdef vector[Vessel] v = vector[Vessel]()
        cdef cyVessel a
        if type(self) is cyVesselNetworkSc:
            if vec is not None and vessels is None:
                self.thisptrDerived = self.thisptr = new VesselNetworkSc(deref(vec.thisptr), p0, rho, Re, dx, Nx,
                                                                         qc, rc)
            if vessels is not None and vec is None:
                self._vessels = vessels
                for a in vessels:
                    v.push_back(deref(a.thisptr))
                self.thisptrDerived = self.thisptr = new VesselNetworkSc(v, p0, rho, Re, dx, Nx, qc, rc)

    def set_boundary_layer_th(self, T, no_cycles):
        self.thisptrDerived.set_boundary_layer_th(T, no_cycles)

    def __dealloc__(self):
        if type(self) is cyVesselNetworkSc:
            del self.thisptrDerived


cpdef double pystd_dev(np.ndarray[np.float64_t, ndim=1] a):
    return std_dev(<double*> a.data, a.size)

cpdef int pytdma(np.ndarray[np.float64_t, ndim=1] a,
           np.ndarray[np.float64_t, ndim=1] b,
           np.ndarray[np.float64_t, ndim=1] c,
           np.ndarray[np.float64_t, ndim=1] d,
           np.ndarray[np.float64_t, ndim=1] out):

    cdef Py_ssize_t siz = d.shape[0]
    return tdma(<double*> a.data, <double*> b.data, <double*> c.data,
                <double*> d.data, <double*> out.data, siz)


cpdef int cpwgrad(np.ndarray[np.float64_t, ndim=1] inarr, np.ndarray[np.float64_t, ndim=1] out,
                  double dx):
    cdef Py_ssize_t siz = inarr.shape[0]
    return grad(<double*> inarr.data, <double*> out.data, dx, siz)


@cython.boundscheck(False)
@cython.wraparound(False)
def cwrite2file(string filename, double[:, ::1] u):
    """
    Function to write results in a file.

    :param filename: the name of the file
    :param u: the solution vetctor
    :return: int
    """
    cdef Py_ssize_t rws
    cdef Py_ssize_t clms
    rws = u.shape[0]
    clms = u.shape[1]
    return write2file(filename, &u[0, 0], rws, clms)

cdef enum CFL:
    CFL_STATUS_OK = 0
    CFL_STATUS_ERROR = -1

cdef enum SOLVER:
    SOLVER_STATUS_OK = 0
    SOLVER_STATUS_ERROR = -1

cdef struct U:
    double** u
    double** u_n
    double** u_star
    double** F_
    double** S_
    double** u_nh_mh
    double** u_nh_ph
    double** F_nh_mh
    double** S_nh_mh
    double** F_nh_ph
    double** F_star
    double** S_star
    double** S_nh_ph
    double dx
    Py_ssize_t x_size

cdef struct Solution:
    double* u_store
    Py_ssize_t size
    Py_ssize_t n

cdef struct cvessel:
    double* r0
    double* f_r0
    double* interpolate_R0p
    double* interpolate_R0m
    double* f_r0_ph
    double* f_r0_mh
    double* df_dr0
    double* df_dx
    double* df_dr0_ph
    double* df_dx_ph
    double* df_dr0_mh
    double* df_dx_mh
    double L, dx, phi
    double* RLC
    Py_ssize_t _size
    Py_ssize_t _rlc_size

cdef struct cParabolicSolObj:
    double** x_hyper
    double* a
    double* b
    double* c
    double* d
    double* a_ph
    double* a_mh
    Py_ssize_t _size


cdef class cPDEm:
    """
    .. math::

        \\frac{\\partial U}{\\partial t} + \\frac{\\partial F}{\\partial x} = S(U)

    :param pylsewave.mesh.VesselNetwork mesh: the vessel network
    :param float rho: the density of fluid medium (e.g. blood, water, etc.)
    :param float mu: blood viscosity
    :param float Re: the Reynolds number
    """
    cdef Pool mem
    cdef object mesh
    cdef double delta, Re, _rho, _mu
    cdef vector[cvessel] vessels
    def __init__(self, mesh, rho, mu, Re=None):
        self.mem = Pool()
        self.mesh = mesh
        if Re is not None:
            self.Re = Re
        if rho is not None:
            self._rho = rho
        if mu is not None:
            self._mu = mu
        cdef Py_ssize_t siz, i
        siz = len(self.mesh.vessels)
        self.vessels.resize(siz)
        # set the vessels
        self.init_vessels(self.vessels)
        self.set_vessels(self.vessels, siz)

    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void pressure_visco(self, double** u, rows, ITYPE_t vessel_index,
     double* out_p):
        raise NotImplementedError("This method is not implemented!")

    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef ITYPE_t CV_f(self, double* a, Py_ssize_t rows, ITYPE_t vessel_index, double* out):
        raise NotImplementedError("This method is not implemented!")

    @cython.cdivision(True)
    cpdef void set_boundary_layer_th(self, double T, Py_ssize_t no_cycles):
        """
        Method to calculate the boundary layer.

        .. math::

            \\delta = \\sqrt{\\frac{\\nu T_{cycle}}{2 \\pi}}

        :param double T: Total period of the simulation
        :param int no_cycles: No of cycles that the simulation will run
        """
        cdef double _nu = self._mu / self._rho
        cdef double T_cycle = T / no_cycles
        self.delta = np.power((_nu * T_cycle) / (2. * np.pi), 0.5)

    @property
    def delta(self):
        return self.delta

    cdef void init_vessels(self, vector[cvessel] &vectorV) except *:
        cdef Py_ssize_t _size = vectorV.size()
        cdef Py_ssize_t _rows, _rlc_rows
        for i in range(_size):
            _rows = self.mesh.vessels[i].r0.shape[0]
            vectorV[i]._size = _rows
            vectorV[i].r0 = <double*>self.mem.alloc(_rows, sizeof(double))
            vectorV[i].f_r0 = <double*>self.mem.alloc(_rows, sizeof(double))
            vectorV[i].interpolate_R0p = <double*>self.mem.alloc(_rows, sizeof(double))
            vectorV[i].interpolate_R0m = <double*>self.mem.alloc(_rows, sizeof(double))
            vectorV[i].f_r0_ph = <double*>self.mem.alloc(_rows, sizeof(double))
            vectorV[i].f_r0_mh = <double*>self.mem.alloc(_rows, sizeof(double))
            vectorV[i].df_dr0 = <double*>self.mem.alloc(_rows, sizeof(double))
            vectorV[i].df_dx = <double*>self.mem.alloc(_rows, sizeof(double))
            vectorV[i].df_dr0_ph = <double*>self.mem.alloc(_rows, sizeof(double))
            vectorV[i].df_dx_ph = <double*>self.mem.alloc(_rows, sizeof(double))
            vectorV[i].df_dr0_mh = <double*>self.mem.alloc(_rows, sizeof(double))
            vectorV[i].df_dx_mh = <double*>self.mem.alloc(_rows, sizeof(double))
            if self.mesh.vessels[i].RLC is not None:
                _rlc_rows = np.array(list(self.mesh.vessels[i].RLC.values())).shape[0]
                vectorV[i]._rlc_size = _rlc_rows
                vectorV[i].RLC = <double*>self.mem.alloc(_rlc_rows, sizeof(double))

    cdef void set_vessels(self, vector[cvessel] &vectorV, Py_ssize_t siz) except *:
        cdef Py_ssize_t i, j, k, rows
        for i in range(siz):
            vectorV[i].L = self.mesh.vessels[i].length
            vectorV[i].dx = self.mesh.vessels[i].dx
            rows = self.mesh.vessels[i].r0.shape[0]
            for j in range(rows):
                vectorV[i].r0[j] = self.mesh.vessels[i].r0[j]
                vectorV[i].f_r0[j] = self.mesh.vessels[i].f_r0[j]
                vectorV[i].interpolate_R0p[j] = self.mesh.vessels[i].interpolate_R0(0.5)[j]
                vectorV[i].interpolate_R0m[j] = self.mesh.vessels[i].interpolate_R0(-0.5)[j]
                vectorV[i].f_r0_ph[j] = self.mesh.vessels[i].f_r0_ph[j]
                vectorV[i].f_r0_mh[j] = self.mesh.vessels[i].f_r0_mh[j]
                vectorV[i].df_dr0[j] = self.mesh.vessels[i].df_dr0[j]
                vectorV[i].df_dx[j] = self.mesh.vessels[i].df_dx[j]
                vectorV[i].df_dr0_ph[j] = self.mesh.vessels[i].df_dr0_ph[j]
                vectorV[i].df_dx_ph[j] = self.mesh.vessels[i].df_dx_ph[j]
                vectorV[i].df_dr0_mh[j] = self.mesh.vessels[i].df_dr0_mh[j]
                vectorV[i].df_dx_mh[j] = self.mesh.vessels[i].df_dx_mh[j]

            if self.mesh.vessels[i].RLC is not None:
                rlc_arr = np.array(list(self.mesh.vessels[i].RLC.values()))
                for k in range(vectorV[i]._rlc_size):
                    vectorV[i].RLC[k] = rlc_arr[k]

    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef ITYPE_t flux(self, double** u, Py_ssize_t rows, ITYPE_t x, ITYPE_t vessel_index, double** out_F) nogil:
        """
        Method to calculate the Flux term.

        :param u: the solution vector u
        :type priority: ndarray[ndim=2, dtype=float]
        :param x: the spatial point index
        :type priority: int
        :param int vessel_index: the unique vessel Id
        :param ndarray[ndim=2, dtype=float] out_F: the flux vector
        :return: int
        """
        cdef double A0, f_r0
        cdef Py_ssize_t siz, i

        siz = rows

        if x == 0:
            for i in range(siz):
                A0 = M_PI*pow(self.vessels[vessel_index].r0[i], 2)
                f_r0 = self.vessels[vessel_index].f_r0[i]
                out_F[0][i] = u[1][i]
                out_F[1][i] = ((u[1][i] * u[1][i]) / u[0][i]) + (f_r0/self._rho)*sqrt(A0*u[0][i])
        elif x == 1:
            for i in range(siz):
                A0 = M_PI*pow(self.vessels[vessel_index].interpolate_R0p[i], 2)
                f_r0 = self.vessels[vessel_index].f_r0_ph[i]
                out_F[0][i] = u[1][i]
                out_F[1][i] = ((u[1][i]*u[1][i]) / u[0][i]) + (f_r0/self._rho)*sqrt(A0*u[0][i])
        elif x == -1:
            for i in range(siz):
                A0 = M_PI*pow(self.vessels[vessel_index].interpolate_R0m[i], 2)
                f_r0 = self.vessels[vessel_index].f_r0_mh[i]
                out_F[0][i] = u[1][i]
                out_F[1][i] = ((u[1][i] * u[1][i]) / u[0][i]) + (f_r0/self._rho)*sqrt(A0*u[0][i])

        return 0

    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cpdef void flux_i(self, double[::1] u, ITYPE_t x, ITYPE_t index, ITYPE_t vessel_index, double[::1] out_f):
        """
        Method to calculate the Flux term.

        :param u: the solution vector u
        :type priority: float or ndarray[ndim=1, dtype=double]
        :param x: the spatial point index
        :type priority: int
        :param int vessel_index: the unique vessel Id
        :param ndarray[ndim=1, dtype=double] out_F: the flux vector
        :return: int
        """
        cdef double A0, f_r0
        cdef Py_ssize_t siz, i

        if x == 0:
            A0 = M_PI*pow(self.vessels[vessel_index].r0[index], 2)
            f_r0 = self.vessels[vessel_index].f_r0[index]
            out_f[0] = u[1]
            out_f[1] = ((u[1]*u[1]) / u[0]) + (f_r0/self._rho)*sqrt(A0*u[0])
        elif x == 1:
            A0 = M_PI*pow(self.vessels[vessel_index].interpolate_R0p[index], 2)
            f_r0 = self.vessels[vessel_index].f_r0_ph[index]
            out_f[0] = u[1]
            out_f[1] = ((u[1]*u[1]) / u[0]) + (f_r0/self._rho)*sqrt(A0*u[0])
        elif x == -1:
            A0 = M_PI*pow(self.vessels[vessel_index].interpolate_R0m[index], 2)
            f_r0 = self.vessels[vessel_index].f_r0_mh[index]
            out_f[0] = u[1]
            out_f[1] = ((u[1]*u[1]) / u[0]) + (f_r0/self._rho)*sqrt(A0*u[0])

    property mesh:
        def __get__(self):
            return self.mesh

        def __set__(self, mesh):
            self.mesh = mesh

    cdef vector[cvessel] get_vessels(self):
        return self.vessels

    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef ITYPE_t source(self, double** u, Py_ssize_t clms, ITYPE_t x, ITYPE_t vessel_index, double** out_S) nogil:
        """
        Method to calculate the Source term.

        :param u: the solution vector u
        :type priority: float or ndarray[ndim=2, dtype=double]
        :param x: the spatial points discretising the 1D segments
        :type priority: float or ndarray[ndim=2, dtype=double]
        :param int vessel_index: the unique vessel Id
        :param ndarray[ndim=2, dtype=double] out_S: the source vector 
        """
        cdef double A0, R, f_r0, df_dr0, df_dx
        cdef Py_ssize_t siz, i
        cdef double _nu = ((self._mu)/self._rho)
        siz = clms

        if x == 0:
            for i in range(siz):
                R = sqrt(u[0][i] / M_PI)
                A0 = M_PI*pow(self.vessels[vessel_index].r0[i], 2)
                f_r0 = self.vessels[vessel_index].f_r0[i]
                df_dr0 = self.vessels[vessel_index].df_dr0[i]
                df_dx = self.vessels[vessel_index].df_dx[i]
                out_S[0][i] = 0.0
                out_S[1][i] = (((-2.*M_PI*R*_nu*u[1][i])/(self.delta*u[0][i])) +
                               (1./self._rho)*(2.*sqrt(u[0][i])*(sqrt(M_PI)*f_r0 + sqrt(A0)*df_dr0) - u[0][i]*df_dr0)*df_dx)
        elif x == 1:
            for i in range(siz):
                R = sqrt(u[0][i] / M_PI)
                A0 = M_PI*pow(self.vessels[vessel_index].interpolate_R0p[i], 2)
                f_r0 = self.vessels[vessel_index].f_r0_ph[i]
                df_dr0 = self.vessels[vessel_index].df_dr0_ph[i]
                df_dx = self.vessels[vessel_index].df_dx_ph[i]
                out_S[0][i] = 0.0
                out_S[1][i] = (((-2.*M_PI*R*_nu*u[1][i]) / (self.delta*u[0][i])) +
                               (1./self._rho)*(2.*sqrt(u[0][i])*(sqrt(M_PI)*f_r0 + sqrt(A0)*df_dr0) - u[0][i] * df_dr0)*df_dx)
        elif x == -1:
            for i in range(siz):
                R = sqrt(u[0][i] / M_PI)
                A0 = M_PI*pow(self.vessels[vessel_index].interpolate_R0m[i], 2)
                f_r0 = self.vessels[vessel_index].f_r0_mh[i]
                df_dr0 = self.vessels[vessel_index].df_dr0_mh[i]
                df_dx = self.vessels[vessel_index].df_dx_mh[i]
                out_S[0][i] = 0.0
                out_S[1][i] = (((-2.*M_PI*_nu*R*u[1][i]) / (self.delta*u[0][i])) +
                               (1./self._rho)*(2.*sqrt(u[0][i])*(sqrt(M_PI)*f_r0 + sqrt(A0)*df_dr0) - u[0][i] * df_dr0)*df_dx)

        return 0

    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cpdef void source_i(self, double[::1] u, ITYPE_t x, ITYPE_t index, ITYPE_t vessel_index, double[::1] out_s):
        """
        Method to calculate the Source term.

        :param u: the solution vector u
        :type priority: float or ndarray[ndim=1, dtype=double]
        :param x: the spatial points discretising the 1D segments
        :type priority: float or ndarray[ndim=1, dtype=double]
        :param int vessel_index: the unique vessel Id
        :param ndarray[ndim=1, dtype=double] out_S: the source vector 
        """
        cdef double A0, R, f_r0, df_dr0, df_dx
        cdef double _nu = ((self._mu)/self._rho)

        if x == 0:
            R = sqrt(u[0] / M_PI)
            A0 = M_PI*pow(self.vessels[vessel_index].r0[index], 2)
            f_r0 = self.vessels[vessel_index].f_r0[index]
            df_dr0 = self.vessels[vessel_index].df_dr0[index]
            df_dx = self.vessels[vessel_index].df_dx[index]
            out_s[0] = 0.0
            out_s[1] = (((-2.*M_PI*_nu*R*u[1]) / (self.delta*u[0])) +
                        (1./self._rho)*(2.*sqrt(u[0])*(sqrt(M_PI)*f_r0 + sqrt(A0)*df_dr0) -
                         u[0] * df_dr0)*df_dx)
        elif x == 1:
            R = sqrt(u[0] / M_PI)
            A0 = M_PI*pow(self.vessels[vessel_index].interpolate_R0p[index], 2)
            f_r0 = self.vessels[vessel_index].f_r0_ph[index]
            df_dr0 = self.vessels[vessel_index].df_dr0_ph[index]
            df_dx = self.vessels[vessel_index].df_dx_ph[index]
            out_s[0] = 0.0
            out_s[1] = (((-2.*M_PI*_nu*R*u[1]) / (self.delta*u[0])) +
                        (1./self._rho)*(2.*sqrt(u[0])*(sqrt(M_PI)*f_r0 + sqrt(A0)*df_dr0) -
                         u[0] * df_dr0)*df_dx)

        elif x == -1:
            R = sqrt(u[0] / M_PI)
            A0 = M_PI*pow(self.vessels[vessel_index].interpolate_R0m[index], 2)
            f_r0 = self.vessels[vessel_index].f_r0_mh[index]
            df_dr0 = self.vessels[vessel_index].df_dr0_mh[index]
            df_dx = self.vessels[vessel_index].df_dx_mh[index]
            out_s[0] = 0.0
            out_s[1] = (((-2.*M_PI*_nu*R*u[1]) / (self.delta*u[0])) +
                        (1./self._rho)*(2.*sqrt(u[0])*(sqrt(M_PI)*f_r0 + sqrt(A0)*df_dr0) -
                         u[0] * df_dr0)*df_dx)

    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef double pressure_i(self, double a, ITYPE_t index, ITYPE_t vessel_index) nogil:
        """
        Method to compute pressure from pressure-area relationship

        .. math::
             p(A) = f(R_0, k) \\left( \\sqrt{\\ 1 - frac{A_0}{A}}\\right)

        :param a: cross-sectional area of the vessel
        :type priority: float or ndarray[ndim=1, type=double]
        :param int index: index of the node in vessel
        :param int vessel_index: the vessel unique Id
        :return: double
        """
        cdef double A0
        cdef double out
        A0 = M_PI*pow(self.vessels[vessel_index].r0[index], 2)
        out = self.vessels[vessel_index].f_r0[index]*(1 - sqrt(A0 / a))
        return out

    @cython.initializedcheck(False)
    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void pressure(self, double[:, ::1] u, ITYPE_t vessel_index, double[::1] out_p):
        """
        Method to compute pressure from pressure-area relationship

        .. math::
             p(A) = f(R_0, k) \\left( \\sqrt{\\ 1 - frac{A_0}{A}}\\right)

        :param u: solution vector
        :type priority: ndarray[ndim=2, type=double]
        :param int index: index of the node in vessel
        :param int vessel_index: the vessel unique Id
        :param ndarray[ndim=1, type=double] out_p: the calculated pressue vector.
        """
        cdef double A0
        cdef Py_ssize_t i, siz
        siz = u.shape[1]
        for i in range(siz):
            A0 = M_PI*pow(self.vessels[vessel_index].r0[i], 2)
            out_p[i] = self.vessels[vessel_index].f_r0[i]*(1 - sqrt(A0 / u[0, i]))

    @cython.initializedcheck(False)
    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef double compute_c_i(self, double a, ITYPE_t index, ITYPE_t vessel_index):
        """
        Method to compute local wave speed

        .. math::
             c(A) = -\\sqrt{\\frac{1}{2 \\rho} f(R_0, k) \\sqrt{\\frac{A_0}{A}}}

        :param A: cross-sectional area of vessel
        :type priority: double
        :param int index: the number or the node to be calculated
        :param vessel_index: the vessel unique Id
        :return: double
        """
        cdef double A0
        cdef double out
        A0 = M_PI*pow(self.vessels[vessel_index].r0[index], 2)
        out = -sqrt((0.5/self._rho)*self.vessels[vessel_index].f_r0[index]*sqrt(A0/a))
        return out

    # @cython.initializedcheck(False)
    # @cython.cdivision(True)
    # @cython.boundscheck(False)
    cpdef void compute_c(self, double[:, ::1] u, ITYPE_t vessel_index, double[::1] out_c):
        """
        Method to compute local wave speed

        .. math::
             c(A) = -\\sqrt{\\frac{1}{2 \\rho} f(R_0, k) \\sqrt{\\frac{A_0}{A}}}

        :param u: solution vector
        :type priority: ndarray[ndim=2, type=double]
        :param int index: the number or the node to be calculated
        :param vessel_index: the vessel unique Id
        :param ndarray[ndim=2, type=double] out_c: the calculated wave speed vector
        """
        cdef double A0
        cdef Py_ssize_t i, siz
        siz = u.shape[1]
        for i in range(siz):
            A0 = M_PI*pow(self.vessels[vessel_index].r0[i], 2)
            out_c[i] = -sqrt((0.5/self._rho)*self.vessels[vessel_index].f_r0[i]*sqrt(A0/u[0, i]))


# ------------------- PDEs from Watanabe and modified by Ragkousis ----------------------- #
cdef class cPDEsWat(cPDEm):
    @cython.initializedcheck(False)
    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef double compute_c_i(self, double a, ITYPE_t index, ITYPE_t vessel_index):
        """
        Method to compute local wave speed

        .. math::
             c(A) = \\sqrt{\\frac{1}{2 \\rho} f(R_0, k) \\sqrt{\\frac{A}{A_0}}}

        :param a: the cross sectional area of the vessel
        :type priority: double
        :param int index: the number or the node to be calculated
        :param vessel_index: the vessel unique Id
        :return: double
        """
        cdef double A0
        cdef double out
        A0 = M_PI*pow(self.vessels[vessel_index].r0[index], 2)
        out = sqrt((0.5/self._rho)*self.vessels[vessel_index].f_r0[index]*sqrt(a/A0))
        return out

    @cython.initializedcheck(False)
    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef void compute_c(self, double[:, ::1] u, ITYPE_t vessel_index, double[::1] out_c):
        """
        Method to compute local wave speed

        .. math::
             c(A) = \\sqrt{\\frac{1}{2 \\rho} f(R_0, k) \\sqrt{\\frac{A}{A_0}}}

        :param u: the solution vector
        :type priority: ndarray[ndim=2, type=double]
        :param int index: the number or the node to be calculated
        :param vessel_index: the vessel unique Id
        :param out_c: the calculated wave speed vector
        :type priority: ndarray[ndim=1, type=double]
        """
        cdef double A0
        cdef Py_ssize_t i, siz
        siz = u.shape[1]
        for i in range(siz):
            A0 = M_PI*pow(self.vessels[vessel_index].r0[i], 2)
            out_c[i] = sqrt((0.5/self._rho)*self.vessels[vessel_index].f_r0[i]*sqrt(u[0, i]/A0))

    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef double pressure_i(self, double a, ITYPE_t index, ITYPE_t vessel_index) nogil:
        """
        Method to compute pressure from pressure-area relationship

        .. math::
             p(A) = f(R_0, k) \\left( \\sqrt{\\ frac{A}{A_0} - 1}\\right)

        :param a: cross sectional area
        :type priority: double
        :param int index: index of the node in vessel
        :param int vessel_index: the vessel unique Id
        :return: double
        """
        cdef double A0
        cdef double out
        A0 = M_PI*pow(self.vessels[vessel_index].r0[index], 2)
        out = self.vessels[vessel_index].f_r0[index]*(sqrt(a / A0) - 1.0)
        return out

    @cython.initializedcheck(False)
    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void pressure(self, double[:, ::1] u, ITYPE_t vessel_index, double[::1] out_p):
        """
        Method to compute pressure from pressure-area relationship

        .. math::
             p(A) = f(R_0, k) \\left( \\sqrt{\\ frac{A}{A_0} - 1}\\right)

        :param u: solution vector
        :type priority: ndarray[ndim=2, type=double]
        :param int index: index of the node in vessel
        :param int vessel_index: the vessel unique Id
        :param ndarray[ndim=1, type=double] out_p: the calculated pressue vector.
        """
        cdef double A0
        cdef Py_ssize_t i, siz
        siz = u.shape[1]
        for i in range(siz):
            A0 = M_PI*pow(self.vessels[vessel_index].r0[i], 2)
            out_p[i] = self.vessels[vessel_index].f_r0[i]*(sqrt(u[0, i] / A0) - 1.0)

    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef ITYPE_t flux(self, double** u, Py_ssize_t rows, ITYPE_t x, ITYPE_t vessel_index, double** out_F) nogil:
        """
        Method to calculate the Flux term.
    
        :param u: the solution vector u
        :type priority: ndarray[ndim=2, dtype=double]
        :param x: the spatial point
        :type priority: int
        :param int vessel_index: the unique vessel Id
        :param ndarray[ndim=2, dtype=double] out_F: the flux vector
        :return: int
        """
        cdef double A0, f_r0
        cdef Py_ssize_t siz, i

        siz = rows

        if x == 0:
            for i in range(siz):
                A0 = M_PI*pow(self.vessels[vessel_index].r0[i], 2)
                f_r0 = self.vessels[vessel_index].f_r0[i]
                out_F[0][i] = u[1][i]
                out_F[1][i] = ((u[1][i]*u[1][i]) / u[0][i]) + (f_r0/(3.0*self._rho))*((u[0][i]**(3/2.))*(A0**(-1/2.)))
        elif x == 1:
            for i in range(siz):
                A0 = M_PI*pow(self.vessels[vessel_index].interpolate_R0p[i], 2)
                f_r0 = self.vessels[vessel_index].f_r0_ph[i]
                out_F[0][i] = u[1][i]
                out_F[1][i] = ((u[1][i]*u[1][i]) / u[0][i]) + (f_r0/(3.0*self._rho))*((u[0][i]**(3/2.))*(A0**(-1/2.)))
        elif x == -1:
            for i in range(siz):
                A0 = M_PI*pow(self.vessels[vessel_index].interpolate_R0m[i], 2)
                f_r0 = self.vessels[vessel_index].f_r0_mh[i]
                out_F[0][i] = u[1][i]
                out_F[1][i] = ((u[1][i]*u[1][i]) / u[0][i]) + (f_r0/(3.0*self._rho))*((u[0][i]**(3/2.))*(A0**(-1/2.)))

        return 0

    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cpdef void flux_i(self, double[::1] u, ITYPE_t x, ITYPE_t index, ITYPE_t vessel_index, double[::1] out_f):
        """
        Method to calculate the Flux term.

        :param u: the solution vector u
        :type priority: ndarray[ndim=1, dtype=double]
        :param x: the spatial point index
        :type priority: int
        :param int vessel_index: the unique vessel Id
        :param ndarray[ndim=1, dtype=double] out_F: the flux vector
        :return: None
        """
        cdef double A0, f_r0
        cdef Py_ssize_t siz, i

        if x == 0:
            A0 = M_PI*pow(self.vessels[vessel_index].r0[index], 2)
            f_r0 = self.vessels[vessel_index].f_r0[index]
            out_f[0] = u[1]
            out_f[1] = ((u[1]*u[1]) / u[0]) + (f_r0/(3.0*self._rho))*((u[0]**(3/2.))*(A0**(-1/2.)))
        elif x == 1:
            A0 = M_PI*pow(self.vessels[vessel_index].interpolate_R0p[index], 2)
            f_r0 = self.vessels[vessel_index].f_r0_ph[index]
            out_f[0] = u[1]
            out_f[1] = ((u[1]*u[1]) / u[0]) + (f_r0/(3.0*self._rho))*((u[0]**(3/2.))*(A0**(-1/2.)))
        elif x == -1:
            A0 = M_PI*pow(self.vessels[vessel_index].interpolate_R0m[index], 2)
            f_r0 = self.vessels[vessel_index].f_r0_mh[index]
            out_f[0] = u[1]
            out_f[1] = ((u[1]*u[1]) / u[0]) + (f_r0/(3.0*self._rho))*((u[0]**(3/2.))*(A0**(-1/2.)))

    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef ITYPE_t source(self, double** u, Py_ssize_t clms, ITYPE_t x, ITYPE_t vessel_index, double** out_S) nogil:
        """
         Method to calculate the Source term.

        :param u: the solution vector u
        :type priority: ndarray[ndim=2, dtype=double]
        :param x: the spatial point index
        :type priority: int
        :param int vessel_index: the unique vessel Id
        :param ndarray[ndim=2, dtype=double] out_S: the source vector 
        """
        cdef double A0, R, f_r0, df_dr0, df_dx
        cdef Py_ssize_t siz, i
        cdef double _nu = ((self._mu)/self._rho)
        siz = clms

        if x == 0:
            for i in range(siz):
                R = sqrt(u[0][i] / M_PI)
                A0 = M_PI*pow(self.vessels[vessel_index].r0[i], 2)
                f_r0 = self.vessels[vessel_index].f_r0[i]
                df_dr0 = self.vessels[vessel_index].df_dr0[i]
                df_dx = self.vessels[vessel_index].df_dx[i]
                out_S[0][i] = 0.0
                out_S[1][i] = (((-2.*M_PI*_nu*R*u[1][i]) / (self.delta*u[0][i])) +
                               (1./self._rho)*(((2*sqrt(M_PI)*f_r0*u[0][i]**(3/2.))/(3.*A0)) - (((2/3.)*(u[0][i]**(3/2.))*A0**(-1/2.)) - u[0][i])*df_dr0)*df_dx)
        elif x == 1:
            for i in range(siz):
                R = sqrt(u[0][i] / M_PI)
                A0 = M_PI*pow(self.vessels[vessel_index].interpolate_R0p[i], 2)
                f_r0 = self.vessels[vessel_index].f_r0_ph[i]
                df_dr0 = self.vessels[vessel_index].df_dr0_ph[i]
                df_dx = self.vessels[vessel_index].df_dx_ph[i]
                out_S[0][i] = 0.0
                out_S[1][i] = (((-2.*M_PI*_nu*R*u[1][i]) / (self.delta*u[0][i])) +
                               (1./self._rho)*(((2*sqrt(M_PI)*f_r0*u[0][i]**(3/2.))/(3.*A0)) - (((2/3.)*(u[0][i]**(3/2.))*A0**(-1/2.)) - u[0][i])*df_dr0)*df_dx)
        elif x == -1:
            for i in range(siz):
                R = sqrt(u[0][i] / M_PI)
                A0 = M_PI*pow(self.vessels[vessel_index].interpolate_R0m[i], 2)
                f_r0 = self.vessels[vessel_index].f_r0_mh[i]
                df_dr0 = self.vessels[vessel_index].df_dr0_mh[i]
                df_dx = self.vessels[vessel_index].df_dx_mh[i]
                out_S[0][i] = 0.0
                out_S[1][i] = (((-2.*M_PI*_nu*R*u[1][i]) / (self.delta*u[0][i])) +
                               (1./self._rho)*(((2*sqrt(M_PI)*f_r0*u[0][i]**(3/2.))/(3.*A0)) - (((2/3.)*(u[0][i]**(3/2.))*A0**(-1/2.)) - u[0][i])*df_dr0)*df_dx)

        return 0

    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cpdef void source_i(self, double[::1] u, ITYPE_t x, ITYPE_t index, ITYPE_t vessel_index, double[::1] out_s):
        """
         Method to calculate the Source term.

        :param u: the solution vector u
        :type priority: ndarray[ndim=1, dtype=double]
        :param x: the spatial point index
        :type priority: int
        :param int vessel_index: the unique vessel Id
        :param ndarray[ndim=1, dtype=double] out_S: the source vector 
        """
        cdef double A0, R, f_r0, df_dr0, df_dx
        cdef double _nu = ((self._mu)/self._rho)

        if x == 0:
            R = sqrt(u[0] / M_PI)
            A0 = M_PI*pow(self.vessels[vessel_index].r0[index], 2)
            f_r0 = self.vessels[vessel_index].f_r0[index]
            df_dr0 = self.vessels[vessel_index].df_dr0[index]
            df_dx = self.vessels[vessel_index].df_dx[index]
            out_s[0] = 0.0
            out_s[1] = (((-2.*M_PI*_nu*R*u[1]) / (self.delta*u[0])) +
                        (1./self._rho)*(((2*sqrt(M_PI)*f_r0*u[0]**(3/2.))/(3.*A0)) - (((2/3.)*(u[0]**(3/2.))*A0**(-1/2.)) - u[0])*df_dr0)*df_dx)
        elif x == 1:
            R = sqrt(u[0] / M_PI)
            A0 = M_PI*pow(self.vessels[vessel_index].interpolate_R0p[index], 2)
            f_r0 = self.vessels[vessel_index].f_r0_ph[index]
            df_dr0 = self.vessels[vessel_index].df_dr0_ph[index]
            df_dx = self.vessels[vessel_index].df_dx_ph[index]
            out_s[0] = 0.0
            out_s[1] = (((-2.*M_PI*_nu*R*u[1]) / (self.delta*u[0])) +
                        (1./self._rho)*(((2*sqrt(M_PI)*f_r0*u[0]**(3/2.))/(3.*A0)) - (((2/3.)*(u[0]**(3/2.))*A0**(-1/2.)) - u[0])*df_dr0)*df_dx)

        elif x == -1:
            R = sqrt(u[0] / M_PI)
            A0 = M_PI*pow(self.vessels[vessel_index].interpolate_R0m[index], 2)
            f_r0 = self.vessels[vessel_index].f_r0_mh[index]
            df_dr0 = self.vessels[vessel_index].df_dr0_mh[index]
            df_dx = self.vessels[vessel_index].df_dx_mh[index]
            out_s[0] = 0.0
            out_s[1] = (((-2.*M_PI*_nu*R*u[1]) / (self.delta*u[0])) +
                        (1./self._rho)*(((2*sqrt(M_PI)*f_r0*u[0]**(3/2.))/(3.*A0)) - (((2/3.)*(u[0]**(3/2.))*A0**(-1/2.)) - u[0])*df_dr0)*df_dx)


# # ------------------- -------------------------------------------------- ----------------------- #
cdef class cPDEsWatVisco(cPDEsWat):
    """
    PDEs class to discretise a system with visco-elastic properties.

    ..note:: see Watanabe at al.
    """
    cdef void set_vessels(self, vector[cvessel] &vectorV, Py_ssize_t siz) except *:
        cdef Py_ssize_t i, j, k, rows
        for i in range(siz):
            vectorV[i].L = self.mesh.vessels[i].length
            vectorV[i].dx = self.mesh.vessels[i].dx
            vectorV[i].phi = self.mesh.vessels[i].phi
            rows = self.mesh.vessels[i].r0.shape[0]
            for j in range(rows):
                vectorV[i].r0[j] = self.mesh.vessels[i].r0[j]
                vectorV[i].f_r0[j] = self.mesh.vessels[i].f_r0[j]
                vectorV[i].interpolate_R0p[j] = self.mesh.vessels[i].interpolate_R0(0.5)[j]
                vectorV[i].interpolate_R0m[j] = self.mesh.vessels[i].interpolate_R0(-0.5)[j]
                vectorV[i].f_r0_ph[j] = self.mesh.vessels[i].f_r0_ph[j]
                vectorV[i].f_r0_mh[j] = self.mesh.vessels[i].f_r0_mh[j]
                vectorV[i].df_dr0[j] = self.mesh.vessels[i].df_dr0[j]
                vectorV[i].df_dx[j] = self.mesh.vessels[i].df_dx[j]
                vectorV[i].df_dr0_ph[j] = self.mesh.vessels[i].df_dr0_ph[j]
                vectorV[i].df_dx_ph[j] = self.mesh.vessels[i].df_dx_ph[j]
                vectorV[i].df_dr0_mh[j] = self.mesh.vessels[i].df_dr0_mh[j]
                vectorV[i].df_dx_mh[j] = self.mesh.vessels[i].df_dx_mh[j]

            if self.mesh.vessels[i].RLC is not None:
                rlc_arr = np.array(list(self.mesh.vessels[i].RLC.values()))
                for k in range(vectorV[i]._rlc_size):
                    vectorV[i].RLC[k] = rlc_arr[k]

    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void pressure_visco(self, double** u, rows, ITYPE_t vessel_index,
     double* out_p):
        """
        Method to compute pressure from pressure-area relationship

        .. math::
             p(A) = f(R_0, k) \\left( \\sqrt{\\frac{A}{A_0} - 1} - C_v \\frac{\\partial q}{\\partial x} \\right)

        :param u: solution vector
        :type priority: ndarray[ndim=2, type=double]
        :param int vessel_index: the vessel unique Id
        :param ndarray[ndim=1, type=double] out_p: the calculated pressue vector.
        """
        cdef Pool mem = Pool()
        cdef double A0
        cdef Py_ssize_t i, siz
        cdef ITYPE_t res
        siz = rows
        cdef double* dqdx = <double*>mem.alloc(siz, sizeof(double))
        cdef double* C_v = <double*>mem.alloc(siz, sizeof(double))
        cdef double* a = <double*>mem.alloc(siz, sizeof(double))
        cdef double* q = <double*>mem.alloc(siz, sizeof(double))
        for i in range(siz):
            a[i] = u[0][i]
            q[i] = u[1][i]
        
        res = grad(q, dqdx, self.vessels[vessel_index].dx, siz)

        res = self.CV_f(a, siz, vessel_index, C_v)
        for i in range(siz):
            A0 = M_PI*pow(self.vessels[vessel_index].r0[i], 2)
            out_p[i] = (self.vessels[vessel_index].f_r0[i]*(sqrt(u[0][i] / A0) - 1.0) -
                        C_v[i]*dqdx[i])

    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef ITYPE_t CV_f(self, double* a, Py_ssize_t rows, ITYPE_t vessel_index, double* out):
        """
        Method to calculate the CV_f term
        
        .. math::
             CV_f = \\frac{2}{3}\\left( \\frac{\\phi w_{th} \\sqrt{\\pi}}{A_0\\sqrt{a}} \\right)
        
        :param a: the cross-sectional area 
        :type priority: ndarray[ndim=1, type=double]
        :param vessel_index: the vessel index
        :type priority: int
        :param out: the calculated vector
        :type priority: ndarray[ndim=1, type=double]
        :return: int
        """
        cdef Pool mem = Pool()
        cdef Py_ssize_t i
        cdef ITYPE_t res
        cdef double A0
        cdef double _phi = self.vessels[vessel_index].phi
        cdef double* w_th = <double*>mem.alloc(rows, sizeof(double))
        res = cPDEsWatVisco.wall_th_x(self.vessels[vessel_index].r0, w_th, rows)
        for i in range(rows):
            A0 = M_PI*pow(self.vessels[vessel_index].r0[i], 2)
            out[i] = (2./3)*((_phi * w_th[i]*sqrt(M_PI))/(A0*sqrt(a[i])))

        return 0

    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @staticmethod
    cdef ITYPE_t wall_th_x(double* r0, double* out, Py_ssize_t rows):
        """
        Method to calculate the wall thickness
        
        :param r0: the un-pressurised vessel radius
        :type priority: ndarray[ndim=1, type=double]
        :param out: the calculated wall thickness
        :type priority: ndarray[ndim=1, type=double]
        :return: int
        """
        cdef Py_ssize_t i
        cdef Py_ssize_t siz = rows
        cdef double alpha = 0.2802
        cdef double beta = -5.053*0.1
        cdef double gamma = 0.1324
        cdef double d = -0.1114*0.1

        for i in range(siz):
            out[i] = r0[i]*(alpha*exp(beta*r0[i]) + gamma*exp(d*r0[i]))

        return 0
# # ------------------- -------------------------------------------------- ----------------------- #

cdef class PDEmCs(cPDEm):
    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef ITYPE_t flux(self, double** u, Py_ssize_t rows, ITYPE_t x, ITYPE_t vessel_index, double** out_F) nogil:
        """
        Method to calculate the Flux term.
    
        :param u: the solution vector u
        :type priority: ndarray[ndim=2, dtype=double]
        :param x: the spatial point
        :type priority: int
        :param int vessel_index: the unique vessel Id
        :param ndarray[ndim=2, dtype=double] out_F: the flux vector
        :return: int
        """
        cdef double A0, f_r0
        cdef Py_ssize_t siz, i

        siz = rows

        if x == 0:
            for i in range(siz):
                A0 = M_PI*pow(self.vessels[vessel_index].r0[i], 2)
                f_r0 = self.vessels[vessel_index].f_r0[i]
                out_F[0][i] = u[1][i]
                out_F[1][i] = ((u[1][i] * u[1][i]) / u[0][i]) + f_r0*sqrt(A0*u[0][i])
        elif x == 1:
            for i in range(siz):
                A0 = M_PI*pow(self.vessels[vessel_index].interpolate_R0p[i], 2)
                f_r0 = self.vessels[vessel_index].f_r0_ph[i]
                out_F[0][i] = u[1][i]
                out_F[1][i] = ((u[1][i]*u[1][i]) / u[0][i]) + f_r0*sqrt(A0*u[0][i])
        elif x == -1:
            for i in range(siz):
                A0 = M_PI*pow(self.vessels[vessel_index].interpolate_R0m[i], 2)
                f_r0 = self.vessels[vessel_index].f_r0_mh[i]
                out_F[0][i] = u[1][i]
                out_F[1][i] = ((u[1][i] * u[1][i]) / u[0][i]) + f_r0* sqrt(A0*u[0][i])

        return 0

    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cpdef void flux_i(self, double[::1] u, ITYPE_t x, ITYPE_t index, ITYPE_t vessel_index, double[::1] out_f):
        """
        Method to calculate the Flux term.
    
        :param u: the solution vector u
        :type priority: ndarray[ndim=1, dtype=double]
        :param x: the spatial point
        :type priority: int
        :param int vessel_index: the unique vessel Id
        :param ndarray[ndim=1, dtype=double] out_F: the flux vector
        :return: int
        """
        cdef double A0, f_r0
        cdef Py_ssize_t siz, i

        if x == 0:
            A0 = M_PI*pow(self.vessels[vessel_index].r0[index], 2)
            f_r0 = self.vessels[vessel_index].f_r0[index]
            out_f[0] = u[1]
            out_f[1] = ((u[1]*u[1]) / u[0]) + f_r0*sqrt(A0*u[0])
        elif x == 1:
            A0 = M_PI*pow(self.vessels[vessel_index].interpolate_R0p[index], 2)
            f_r0 = self.vessels[vessel_index].f_r0_ph[index]
            out_f[0] = u[1]
            out_f[1] = ((u[1]*u[1]) / u[0]) + f_r0*sqrt(A0*u[0])
        elif x == -1:
            A0 = M_PI*pow(self.vessels[vessel_index].interpolate_R0m[index], 2)
            f_r0 = self.vessels[vessel_index].f_r0_mh[index]
            out_f[0] = u[1]
            out_f[1] = ((u[1]*u[1]) / u[0]) + f_r0*sqrt(A0*u[0])

    property mesh:
        def __get__(self):
            return self.mesh

    cdef vector[cvessel] get_vessels(self):
        return self.vessels

    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef ITYPE_t source(self, double** u, Py_ssize_t clms, ITYPE_t x, ITYPE_t vessel_index, double** out_S) nogil:
        """
         Method to calculate the Source term.

        :param u: the solution vector u
        :type priority: ndarray[ndim=2, dtype=double]
        :param x: the spatial point index
        :type priority: int
        :param int vessel_index: the unique vessel Id
        :param ndarray[ndim=2, dtype=double] out_S: the source vector 
        """
        cdef double A0, R, f_r0, df_dr0, df_dx
        cdef Py_ssize_t siz, i

        siz = clms

        if x == 0:
            for i in range(siz):
                R = sqrt(u[0][i] / M_PI)
                A0 = M_PI*pow(self.vessels[vessel_index].r0[i], 2)
                f_r0 = self.vessels[vessel_index].f_r0[i]
                df_dr0 = self.vessels[vessel_index].df_dr0[i]
                df_dx = self.vessels[vessel_index].df_dx[i]
                out_S[0][i] = 0.0
                out_S[1][i] = (((-2.*M_PI*R*u[1][i])/(self.delta*self.Re*u[0][i])) +
                               (2.*sqrt(u[0][i])*(sqrt(M_PI)*f_r0 + sqrt(A0)*df_dr0) - u[0][i]*df_dr0)*df_dx)
        elif x == 1:
            for i in range(siz):
                R = sqrt(u[0][i] / M_PI)
                A0 = M_PI*pow(self.vessels[vessel_index].interpolate_R0p[i], 2)
                f_r0 = self.vessels[vessel_index].f_r0_ph[i]
                df_dr0 = self.vessels[vessel_index].df_dr0_ph[i]
                df_dx = self.vessels[vessel_index].df_dx_ph[i]
                out_S[0][i] = 0.0
                out_S[1][i] = (((-2.*M_PI*R*u[1][i]) / (self.delta*self.Re*u[0][i])) +
                               (2.*sqrt(u[0][i])*(sqrt(M_PI)*f_r0 + sqrt(A0)*df_dr0) - u[0][i] * df_dr0)*df_dx)
        elif x == -1:
            for i in range(siz):
                R = sqrt(u[0][i] / M_PI)
                A0 = M_PI*pow(self.vessels[vessel_index].interpolate_R0m[i], 2)
                f_r0 = self.vessels[vessel_index].f_r0_mh[i]
                df_dr0 = self.vessels[vessel_index].df_dr0_mh[i]
                df_dx = self.vessels[vessel_index].df_dx_mh[i]
                out_S[0][i] = 0.0
                out_S[1][i] = (((-2.*M_PI*R*u[1][i]) / (self.delta*self.Re*u[0][i])) +
                               (2.*sqrt(u[0][i])*(sqrt(M_PI)*f_r0 + sqrt(A0)*df_dr0) - u[0][i] * df_dr0)*df_dx)

        return 0

    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cpdef void source_i(self, double[::1] u, ITYPE_t x, ITYPE_t index, ITYPE_t vessel_index, double[::1] out_s):
        """
         Method to calculate the Source term.

        :param u: the solution vector u
        :type priority: ndarray[ndim=1, dtype=double]
        :param x: the spatial point index
        :type priority: int
        :param int vessel_index: the unique vessel Id
        :param ndarray[ndim=1, dtype=double] out_S: the source vector 
        """
        cdef double A0, R, f_r0, df_dr0, df_dx

        if x == 0:
            R = sqrt(u[0] / M_PI)
            A0 = M_PI*pow(self.vessels[vessel_index].r0[index], 2)
            f_r0 = self.vessels[vessel_index].f_r0[index]
            df_dr0 = self.vessels[vessel_index].df_dr0[index]
            df_dx = self.vessels[vessel_index].df_dx[index]
            out_s[0] = 0.0
            out_s[1] = (((-2.*M_PI*R*u[1]) / (self.delta*self.Re*u[0])) +
                        (2.*sqrt(u[0])*(sqrt(M_PI)*f_r0 + sqrt(A0)*df_dr0) -
                         u[0] * df_dr0)*df_dx)
        elif x == 1:
            R = sqrt(u[0] / M_PI)
            A0 = M_PI*pow(self.vessels[vessel_index].interpolate_R0p[index], 2)
            f_r0 = self.vessels[vessel_index].f_r0_ph[index]
            df_dr0 = self.vessels[vessel_index].df_dr0_ph[index]
            df_dx = self.vessels[vessel_index].df_dx_ph[index]
            out_s[0] = 0.0
            out_s[1] = (((-2.*M_PI*R*u[1]) / (self.delta*self.Re*u[0])) +
                        (2.*sqrt(u[0])*(sqrt(M_PI)*f_r0 + sqrt(A0)*df_dr0) -
                         u[0] * df_dr0)*df_dx)

        elif x == -1:
            R = sqrt(u[0] / M_PI)
            A0 = M_PI*pow(self.vessels[vessel_index].interpolate_R0m[index], 2)
            f_r0 = self.vessels[vessel_index].f_r0_mh[index]
            df_dr0 = self.vessels[vessel_index].df_dr0_mh[index]
            df_dx = self.vessels[vessel_index].df_dx_mh[index]
            out_s[0] = 0.0
            out_s[1] = (((-2.*M_PI*R*u[1]) / (self.delta*self.Re*u[0])) +
                        (2.*sqrt(u[0])*(sqrt(M_PI)*f_r0 + sqrt(A0)*df_dr0) -
                         u[0] * df_dr0)*df_dx)


    @cython.initializedcheck(False)
    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef double compute_c_i(self, double a, ITYPE_t index, ITYPE_t vessel_index):
        """
        Method to compute local wave speed

        .. math::
             c(A) = -\\sqrt{\\frac{1}{2 \\rho} f(R_0, k) \\sqrt{\\frac{A_0)}{A}}}

        :param a: the cross-sectional area of the vessel
        :type priority: double
        :param int index: the number or the node to be calculated
        :param vessel_index: the vessel unique Id
        :return: double
        """
        cdef double A0
        cdef double out
        A0 = M_PI*pow(self.vessels[vessel_index].r0[index], 2)
        out = -sqrt(0.5*self.vessels[vessel_index].f_r0[index]*sqrt(A0/a))
        return out

    # @cython.initializedcheck(False)
    # @cython.cdivision(True)
    # @cython.boundscheck(False)
    cpdef void compute_c(self, double[:, ::1] u, ITYPE_t vessel_index, double[::1] out_c):
        """
        Method to compute local wave speed

        .. math::
             c(A) = -\\sqrt{\\frac{1}{2 \\rho} f(R_0, k) \\sqrt{\\frac{A_0)}{A}}}

        :param u: the solution vector
        :type priority: ndarray[ndim=2, type=double]
        :param int index: the number or the node to be calculated
        :param vessel_index: the vessel unique Id
        :param out_c: the calculated wave speed vector
        :type priority: ndarray[ndim=1, type=double]
        """
        cdef double A0
        cdef Py_ssize_t i, siz
        siz = u.shape[1]
        for i in range(siz):
            A0 = M_PI*pow(self.vessels[vessel_index].r0[i], 2)
            out_c[i] = -sqrt(0.5*self.vessels[vessel_index].f_r0[i]*sqrt(A0/u[0, i]))


cdef class cBCs:
    """
    Class to set boundary conditions. Inlet, outlet, bifs and conjs.

    :param pylsewave.cynum.cPDEm inpdes: the discretised system of pdes.
    :param inletfunc: a python function for inlet BCs.
    """
    cdef cPDEm pdes
    cdef object mesh
    cdef object inlet_fun
    cdef vector[cvessel] vessels
    cdef double _rho
    def __init__(self, cPDEm inpdes, inletfunc):
        self.mesh = inpdes.mesh
        self._rho = inpdes._rho
        self.pdes = inpdes
        self.inlet_fun = inletfunc
        self.vessels = inpdes.get_vessels()

    cdef void set_vessels(self, vector[cvessel] invessels):
        self.vessels = invessels

    cdef cPDEm get_pdes(self):
        return self.pdes

    cdef void set_pdes(self, cPDEm pdes):
        self.pdes = pdes

    # @cython.cdivision(True)
    # @cython.boundscheck(False)
    cdef void U_0(self, double** u, Py_ssize_t clms, double t, double dx, double dt, ITYPE_t vessel_index, double** out):
 #       raise NotImplementedError("Method not implemented")
        cdef Pool mem = Pool()
        cdef double theta, dt2, _q0_nh, _q0_n, q12_nh, A0
        cdef double** S_ = <double**>mem.alloc(2, sizeof(double*))
        cdef double** F_ = <double**>mem.alloc(2, sizeof(double*))
        cdef double[::1] _A, _q, u12_nh
        cdef Py_ssize_t rws, i, j, k
        rws = 2
        
        for k in range(2):
            F_[k] = <double*>mem.alloc(clms, sizeof(double))
            S_[k] = <double*>mem.alloc(clms, sizeof(double))

        _A = np.zeros(clms, dtype=np.float)
        _q = np.zeros(clms, dtype=np.float)
        u12_nh = np.zeros(2, dtype=np.float)

        for j in range(clms):
            _A[j] = u[0][j]
            _q[j] = u[1][j]

        theta = dt/dx
        dt2 = dt/2.

        _q0_nh = self.inlet_fun(t - dt / 2.)
        _q0_n = self.inlet_fun(t)


        self.pdes.flux(u, clms, 0, vessel_index, F_)
        self.pdes.source(u, clms, 0, vessel_index, S_)

        # U[1/2, 1/2]
        for i in range(rws):
            u12_nh[i] = ((u[i][1] + u[i][0]) / 2. -
                         0.5 * theta * (F_[i][1] - F_[i][0]) +
                         0.5 * dt2 * (S_[i][1] + S_[i][0]))
        # q[1/2, 1/2]
        q12_nh = u12_nh[1]
        # A[n+1, 0]
        A0 = _A[0] - (2. * theta) * (q12_nh - _q0_nh)

        out[0][0] = A0
        out[1][0] = _q0_n

    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef ITYPE_t U_L(self, double** u, Py_ssize_t clms, double t, double dx, double dt,
                  ITYPE_t vessel_index, double** out) nogil except *:
        cdef double theta, dt2, out_a, out_q, p0, p_out, x
        cdef Py_ssize_t rws, i, j
        cdef double q_m1
        rws = 2
        theta = dt/dx
        dt2 = dt/2.
        q_m1 = out[1][clms-2]

        p_out = self.pdes.pressure_i(u[0][clms-1], clms-1, vessel_index)
        p0 = self.pdes.pressure_i(u[0][clms-1], clms-1, vessel_index)

        cdef Py_ssize_t k = 0
        cdef double R1, Rt, Ct, p_old

        R1 = self.vessels[vessel_index].RLC[0]
        Rt = self.vessels[vessel_index].RLC[1]
        Ct = self.vessels[vessel_index].RLC[2]

        x = (dt / (R1*Rt*Ct))
        while k < 1000:
            p_old = p0
            out_q = (x*p_out - x*(R1 + Rt)*u[1][clms-1] +
                     (p0 - p_out) / R1 + u[1][clms-1])

            out_a = u[0][clms-1] - theta*(out_q - q_m1)
            p0 = self.pdes.pressure_i(out_a, index=clms-1, vessel_index=vessel_index)
            if fabs(p_old - p0) < 1e-7:
                break
            k += 1

        out[0][clms-1] = out_a
        out[1][clms-1] = out_q

        return STATUS_OK

    # @cython.boundscheck(False)
    cpdef void I(self, double x, ITYPE_t index, ITYPE_t vessel_index, double[::1] out):
        cdef double A0
        A0 = M_PI*pow(self.vessels[vessel_index].r0[index], 2)
        out[0] = A0
        out[1] = 0.0

    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef ITYPE_t bifurcation_R(self, double[::1] x, double[::1] u,
                           double dt, ITYPE_t[::1] vessel_index_list, double[::1] out):
        cdef double x1, x2, x3, x4, x5, x6
        cdef ITYPE_t p, d1, d2
        cdef double A_p, q_p, A_d1, q_d1, A_d2, q_d2
        cdef double rho_
        cdef Py_ssize_t last_i, first_i
        cdef double w1, w2, x_i, W

        x1, x2, x3, x4, x5, x6 = x[0], x[1], x[2], x[3], x[4], x[5]
        p, d1, d2 = vessel_index_list[0], vessel_index_list[1], vessel_index_list[2]
        A_p, q_p, A_d1, q_d1, A_d2, q_d2 = u[0], u[1], u[2], u[3], u[4], u[5]
        rho_ = self._rho

        last_i = self.vessels[p]._size - 1
        first_i = 0

        out[0] = -x2 + x4 + x6
        out[1] = (-0.5*rho_*(x2/x1)**2 - self.vessels[p].f_r0[last_i]*(1 -
                  sqrt(M_PI*self.vessels[p].r0[last_i]**2 / x1)) + 0.5*rho_*(x4/x3)**2 +
                  self.vessels[d1].f_r0[first_i]*(1 - sqrt(M_PI*self.vessels[d1].r0[0]**2 / x3)))
        out[2] = (-0.5*rho_*(x2/x1)**2 - self.vessels[p].f_r0[last_i]*(1 -
                  sqrt(M_PI*self.vessels[p].r0[last_i]**2 / x1)) + 0.5*rho_*(x6/x5)**2 +
                  self.vessels[d2].f_r0[0]*(1 - sqrt(M_PI*self.vessels[d2].r0[0]**2 / x5)))
        out[3] = (-x2/x1 - 4*self.pdes.compute_c_i(x1, index=last_i, vessel_index=p) +
                  (q_p/A_p + 4*self.pdes.compute_c_i(A_p, index=last_i-1, vessel_index=p)))
        out[4] = (-x4/x3 + 4*self.pdes.compute_c_i(x3, index=0, vessel_index=d1) +
                  (q_d1/A_d1 - 4*self.pdes.compute_c_i(A_d1, index=1, vessel_index=d1)))
        out[5] = (-x6/x5 + 4*self.pdes.compute_c_i(x5, index=0, vessel_index=d2) +
                  (q_d2/A_d2 - 4*self.pdes.compute_c_i(A_d2, index=1, vessel_index=d2)))

        return STATUS_OK

    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef ITYPE_t bifurcation_Jr(self, double[::1] x, ITYPE_t[::1] vessel_index_list,
                            double[:, ::1] out):
        cdef double x1, x2, x3, x4, x5, x6, rho_
        cdef ITYPE_t p, d1, d2
        cdef Py_ssize_t last_i, first_i

        x1, x2, x3, x4, x5, x6 = x[0], x[1], x[2], x[3], x[4], x[5]
        p, d1, d2 = vessel_index_list[0], vessel_index_list[1], vessel_index_list[2]

        last_i = self.vessels[p]._size - 1
        first_i = 0

        rho_ = self._rho
        out[0, 0] = 0.
        out[0, 1] = -1.
        out[0, 2] = 0.
        out[0, 3] = 1.
        out[0, 4] = 0.
        out[0, 5] = 1.
        out[1, 0] = rho_*(x2**2 / x1**3) - self.vessels[p].f_r0[last_i]*(M_PI*self.vessels[p].r0[last_i]**2)**0.5 / (2*x1**(3/2.))
        out[1, 1] = -rho_*x2 / x1**2
        out[1, 2] = -rho_*(x4**2 / x3**3) + self.vessels[d1].f_r0[0]*(M_PI*self.vessels[d1].r0[0]**2)**0.5 / (2*x3**(3/2.))
        out[1, 3] = rho_*x4 / x3**2
        out[1, 4] = 0.
        out[1, 5] = 0.
        out[2, 0] = rho_*(x2**2 / x1**3) - self.vessels[p].f_r0[last_i]*(M_PI*self.vessels[p].r0[last_i]**2)**0.5 / (2*x1**(3/2.))
        out[2, 1] = -rho_*x2 / x1**2
        out[2, 2] = 0.
        out[2, 3] = 0.
        out[2, 4] = -rho_*(x6**2 / x5**3) + self.vessels[d2].f_r0[0]*(M_PI*self.vessels[d2].r0[0]**2)**0.5 / (2*x5**(3/2.))
        out[2, 5] = rho_*x6 / x5**2
        out[3, 0] = (x2 / x1**2) - 4*(1/4.*sqrt((0.5/rho_)*self.vessels[p].f_r0[last_i]*sqrt(M_PI*self.vessels[p].r0[last_i]**2)))*x1**(-5/4.)
        out[3, 1] = -1./x1
        out[3, 2] = 0.
        out[3, 3] = 0.0
        out[3, 4] = 0.
        out[3, 5] = 0.
        out[4, 0] = 0.0
        out[4, 1] = 0.
        out[4, 2] = (x4 / x3**2) + 4*(1/4.*sqrt((0.5/rho_)*self.vessels[d1].f_r0[0]*sqrt(M_PI*self.vessels[d1].r0[0]**2)))*x3**(-5/4.)
        out[4, 3] = -1./x3
        out[4, 4] = 0.
        out[4, 5] = 0.
        out[5, 0] = 0.
        out[5, 1] = 0.
        out[5, 2] = 0.
        out[5, 3] = 0.
        out[5, 4] = (x6 / x5**2) + 4*(1/4.*sqrt((0.5/rho_)*self.vessels[d2].f_r0[0]*sqrt(M_PI*self.vessels[d2].r0[0]**2)))*x5**(-5/4.)
        out[5, 5] = -1./x5

        return STATUS_OK

    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cpdef ITYPE_t conjuction_R(self, double[::1] x, double[::1] u, double dt, ITYPE_t[::1] vessel_index_list, double[::1] out):
        cdef double x1, x2, x3, x4, rho_
        cdef double A_d1, q_d1, A_d2, q_d2
        cdef ITYPE_t d1, d2
        cdef Py_ssize_t last_i, first_i
        cdef double w1, w2, x_i, W

        x1, x2, x3, x4 = x[0], x[1], x[2], x[3]
        d1, d2 = vessel_index_list[0], vessel_index_list[1]

        #this is better to be improved and instead of shape - 1, last index should be shape
        last_i = self.vessels[d1]._size - 1
        first_i = 0

        A_d1, q_d1, A_d2, q_d2 = u[0], u[1], u[2], u[3]
        rho_ = self._rho

        out[0] = -x2 + x4
        out[1] = (-0.5*rho_*(x2/x1)**2 - self.vessels[d1].f_r0[last_i]*(1 -
                  sqrt(M_PI*self.vessels[d1].r0[last_i]**2 / x1)) + 0.5*rho_*(x4/x3)**2 +
                  self.vessels[d2].f_r0[0]*(1 - sqrt(M_PI*self.vessels[d2].r0[0]**2 / x3)))
        out[2] = (-x2/x1 - 4*self.pdes.compute_c_i(x1, index=last_i, vessel_index=d1) +
                  (q_d1/A_d1 + 4*self.pdes.compute_c_i(A_d1, index=last_i-1, vessel_index=d1)))
        out[3] = (-x4/x3 + 4*self.pdes.compute_c_i(x3, index=0, vessel_index=d2) +
                  (q_d2/A_d2 - 4*self.pdes.compute_c_i(A_d2, index=1, vessel_index=d2)))


        return STATUS_OK

    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cpdef ITYPE_t conjuction_Jr(self, double[::1] x, ITYPE_t[::1] vessel_index_list, double[:,::1] out):
        cdef double x1, x2, x3, x4, rho_
        cdef ITYPE_t d1, d2
        cdef Py_ssize_t last_i, first_i

        x1, x2, x3, x4 = x[0], x[1], x[2], x[3]
        d1, d2 = vessel_index_list[0], vessel_index_list[1]
        rho_ = self._rho
        last_i = self.vessels[d1]._size - 1

        out[0, 0] = 0.
        out[0, 1] = -1.
        out[0, 2] = 0.
        out[0, 3] = 1.
        out[1, 0] = rho_*(x2**2 / x1**3) - self.vessels[d1].f_r0[last_i]*(M_PI*self.vessels[d1].r0[last_i]**2)**0.5 / (2*x1**(3/2.))
        out[1, 1] = -rho_*x2 / x1**2
        out[1, 2] = -rho_*(x4**2 / x3**3) + self.vessels[d2].f_r0[0]*(M_PI*self.vessels[d2].r0[0]**2)**0.5 / (2*x3**(3/2.))
        out[1, 3] = rho_*x4 / x3**2
        out[2, 0] = (x2 / x1**2) - 4*(1/4.*sqrt((0.5/rho_)*self.vessels[d1].f_r0[last_i]*sqrt(M_PI*self.vessels[d1].r0[last_i]**2)))*x1**(-5/4.)
        out[2, 1] = -1./x1
        out[2, 2] = 0.
        out[2, 3] = 0.
        out[3, 0] = 0.
        out[3, 1] = 0.
        out[3, 2] = (x4 / x3**2) + 4*(1/4.*sqrt((0.5/rho_)*self.vessels[d2].f_r0[0]*sqrt(M_PI*self.vessels[d2].r0[0]**2)))*x3**(-5/4.)
        out[3, 3] = -1./x3

        return STATUS_OK

# # ------------------- \begin{BCs from Watanabe and modified by Ragkousis} ----------------------- #
cdef class cBCsWat(cBCs):
    cdef void U_0(self, double** u, Py_ssize_t clms, double t, double dx, double dt, ITYPE_t vessel_index, double** out):
        theta = dt / dx
        dt2 = dt / 2.
        p_pres = self.inlet_fun(t)

        A0 = np.pi*self.mesh.vessels[vessel_index].r0[0]**2.0
        f = self.mesh.vessels[vessel_index].f_r0[0]

        A = A0*((p_pres /f) + 1)**2.0

        W2_0 = u[1][0]/u[0][0] - 4*self.pdes.compute_c_i(u[0][0], index=0, vessel_index=vessel_index)
        W2_1 = u[1][1]/u[0][1] - 4*self.pdes.compute_c_i(u[0][1], index=1, vessel_index=vessel_index)
        l2 = u[1][0]/u[0][0] - self.pdes.compute_c_i(u[0][0], index=0, vessel_index=vessel_index)
        x = -l2*dt

        W2 = linear_extrapolation(x, 0.0, dx, W2_0, W2_1)
        q = A*(W2 + 4*self.pdes.compute_c_i(A, index=0, vessel_index=vessel_index))

        out[0][0] = A
        out[1][0] = q

    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef ITYPE_t bifurcation_R(self, double[::1] x, double[::1] u, double dt, ITYPE_t[::1] vessel_index_list, double[::1] out):
        cdef double x1, x2, x3, x4, x5, x6
        cdef ITYPE_t p, d1, d2
        cdef double A_p, q_p, A_d1, q_d1, A_d2, q_d2
        cdef double rho_
        cdef Py_ssize_t last_i, first_i
        cdef double w1, w2, x_i, W

        x1, x2, x3, x4, x5, x6 = x[0], x[1], x[2], x[3], x[4], x[5]
        p, d1, d2 = vessel_index_list[0], vessel_index_list[1], vessel_index_list[2]
        A_p, q_p, A_d1, q_d1, A_d2, q_d2 = u[0], u[1], u[2], u[3], u[4], u[5]

        last_i = self.vessels[p]._size - 1
        first_i = 0

        out[0] = -x2 + x4 + x6
        out[1] = (-0.5*self._rho*(x2/x1)**2 - self.vessels[p].f_r0[last_i]*(
                  sqrt(x1 / (M_PI*self.vessels[p].r0[last_i]**2)) - 1) + 0.5*self._rho*(x4/x3)**2 +
                  self.vessels[d1].f_r0[0]*(sqrt(x3 / (M_PI*self.vessels[d1].r0[0]**2)) - 1))
        out[2] = (-0.5*self._rho*(x2/x1)**2 - self.vessels[p].f_r0[last_i]*(
                  sqrt(x1 / (M_PI*self.vessels[p].r0[last_i]**2)) - 1) + 0.5*self._rho*(x6/x5)**2 +
                  self.vessels[d2].f_r0[0]*(sqrt(x5 / (M_PI*self.vessels[d2].r0[0]**2)) - 1))
        out[3] = (-x2/x1 - 4*self.pdes.compute_c_i(x1, index=last_i, vessel_index=p) +
                  (q_p/A_p + 4*self.pdes.compute_c_i(A_p, index=last_i-1, vessel_index=p)))
        out[4] = (-x4/x3 + 4*self.pdes.compute_c_i(x3, index=0, vessel_index=d1) +
                  (q_d1/A_d1 - 4*self.pdes.compute_c_i(A_d1, index=1, vessel_index=d1)))
        out[5] = (-x6/x5 + 4*self.pdes.compute_c_i(x5, index=0, vessel_index=d2) +
                  (q_d2/A_d2 - 4*self.pdes.compute_c_i(A_d2, index=1, vessel_index=d2)))

        return STATUS_OK

    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef ITYPE_t bifurcation_Jr(self, double[::1] x, ITYPE_t[::1] vessel_index_list, double[:, ::1] out):
        cdef double x1, x2, x3, x4, x5, x6, rho_
        cdef ITYPE_t p, d1, d2
        cdef Py_ssize_t last_i, first_i

        x1, x2, x3, x4, x5, x6 = x[0], x[1], x[2], x[3], x[4], x[5]
        p, d1, d2 = vessel_index_list[0], vessel_index_list[1], vessel_index_list[2]

        last_i = self.vessels[p]._size - 1
        first_i = 0

        out[0, 0] = 0.
        out[0, 1] = -1.
        out[0, 2] = 0.
        out[0, 3] = 1.
        out[0, 4] = 0.
        out[0, 5] = 1.
        out[1, 0] = self._rho*(x2**2 / x1**3) - 0.5*self.vessels[p].f_r0[last_i]*((M_PI*self.vessels[p].r0[last_i]**2)**(-0.5))*(x1**(-1/2.))
        out[1, 1] =  -self._rho*x2 / x1**2
        out[1, 2] = -self._rho*(x4**2 / x3**3) + 0.5*self.vessels[d1].f_r0[0]*((M_PI*self.vessels[d1].r0[0]**2)**(-0.5))*(x3**(-1/2.))
        out[1, 3] = self._rho*x4 / x3**2
        out[1, 4] = 0.
        out[1, 5] = 0.
        out[2, 0] = self._rho*(x2**2 / x1**3) - 0.5*self.vessels[p].f_r0[last_i]*((M_PI*self.vessels[p].r0[last_i]**2)**(-0.5))*(x1**(-1/2.))
        out[2, 1] = -self._rho*x2 / x1**2
        out[2, 2] = 0.
        out[2, 3] = 0.
        out[2, 4] = -self._rho*(x6**2 / x5**3) + 0.5*self.vessels[d2].f_r0[0]*((M_PI*self.vessels[d2].r0[0]**2)**(-0.5))*(x5**(-1/2.))
        out[2, 5] = self._rho*x6 / x5**2
        out[3, 0] = (x2 / x1**2) - 4*(1/4.*sqrt((0.5/self._rho)*self.vessels[p].f_r0[last_i] / sqrt(M_PI*self.vessels[p].r0[last_i]**2)))*x1**(-3/4.)
        out[3, 1] = -1./x1
        out[3, 2] = 0.0
        out[3, 3] = 0.0
        out[3, 4] = 0.0
        out[3, 5] = 0.0
        out[4, 0] = 0.0
        out[4, 1] = 0.0
        out[4, 2] = (x4 / x3**2) + 4*(1/4.*sqrt((0.5/self._rho)*self.vessels[d1].f_r0[0] / sqrt(M_PI*self.vessels[d1].r0[0]**2)))*x3**(-3/4.)
        out[4, 3] = -1./x3
        out[4, 4] = 0.
        out[4, 5] = 0.
        out[5, 0] = 0.
        out[5, 1] = 0.
        out[5, 2] = 0.
        out[5, 3] = 0.
        out[5, 4] = (x6 / x5**2) + 4*(1/4.*sqrt((0.5/self._rho)*self.vessels[d2].f_r0[0] / sqrt(M_PI*self.vessels[d2].r0[0]**2)))*x5**(-3/4.)
        out[5, 5] = -1./x5

        return STATUS_OK

    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cpdef ITYPE_t conjuction_R(self, double[::1] x, double[::1] u, double dt, ITYPE_t[::1] vessel_index_list, double[::1] out):
        cdef double x1, x2, x3, x4, rho_
        cdef double A_d1, q_d1, A_d2, q_d2
        cdef ITYPE_t d1, d2
        cdef Py_ssize_t last_i, first_i
        cdef double w1, w2, x_i, W

        x1, x2, x3, x4 = x[0], x[1], x[2], x[3]
        d1, d2 = vessel_index_list[0], vessel_index_list[1]

        #this is better to be improved and instead of shape - 1, last index should be shape
        last_i = self.vessels[d1]._size - 1
        first_i = 0

        A_d1, q_d1, A_d2, q_d2 = u[0], u[1], u[2], u[3]

        out[0] = -x2 + x4
        out[1] = (-0.5*self._rho*(x2/x1)**2 - self.vessels[d1].f_r0[last_i]*(
                  sqrt(x1 / (M_PI*self.vessels[d1].r0[last_i]**2)) - 1) + 0.5*self._rho*(x4/x3)**2 +
                  self.vessels[d2].f_r0[0]*(sqrt(x3/ (M_PI*self.vessels[d2].r0[0]**2)) - 1))
        out[2] = (-x2/x1 - 4*self.pdes.compute_c_i(x1, index=last_i, vessel_index=d1) +
                  (q_d1/A_d1 + 4*self.pdes.compute_c_i(A_d1, index=last_i-1, vessel_index=d1)))
        out[3] = (-x4/x3 + 4*self.pdes.compute_c_i(x3, index=0, vessel_index=d2) +
                  (q_d2/A_d2 - 4*self.pdes.compute_c_i(A_d2, index=1, vessel_index=d2)))

        return STATUS_OK

    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cpdef ITYPE_t conjuction_Jr(self, double[::1] x, ITYPE_t[::1] vessel_index_list, double[:,::1] out):
        cdef double x1, x2, x3, x4, rho_
        cdef ITYPE_t d1, d2
        cdef Py_ssize_t last_i, first_i

        x1, x2, x3, x4 = x[0], x[1], x[2], x[3]
        d1, d2 = vessel_index_list[0], vessel_index_list[1]

        last_i = self.vessels[d1]._size - 1

        out[0, 0] = 0.
        out[0, 1] = -1.
        out[0, 2] = 0.
        out[0, 3] = 1.
        out[1, 0] = self._rho*(x2**2 / x1**3) - 0.5*self.vessels[d1].f_r0[last_i]*((M_PI*self.vessels[d1].r0[last_i]**2)**(-0.5))*(x1**(-1/2.))
        out[1, 1] = -self._rho*x2 / x1**2
        out[1, 2] = -self._rho*(x4**2 / x3**3) + 0.5*self.vessels[d2].f_r0[0]*((M_PI*self.vessels[d2].r0[0]**2)**(-0.5))*(x3**(-1/2.))
        out[1, 3] =  self._rho*x4 / x3**2
        out[2, 0] = (x2 / x1**2) - 4*(1/4.*sqrt((0.5/self._rho)*self.vessels[d1].f_r0[last_i] / sqrt(M_PI*self.vessels[d1].r0[last_i]**2)))*x1**(-3/4.)
        out[2, 1] = -1./x1
        out[2, 2] = 0.
        out[2, 3] = 0.
        out[3, 0] = 0.
        out[3, 1] = 0.
        out[3, 2] = (x4 / x3**2) + 4*(1/4.*sqrt((0.5/self._rho)*self.vessels[d2].f_r0[0] / sqrt(M_PI*self.vessels[d2].r0[0]**2)))*x3**(-3/4.)
        out[3, 3] = -1./x3

        return STATUS_OK


# # ------------------- \end{BCs from Watanabe and modified by Ragkousis} ------------------------- #

cdef class BCsSc(cBCs):
    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef ITYPE_t bifurcation_R(self, double[::1] x, double[::1] u, double dt, ITYPE_t[::1] vessel_index_list, double[::1] out):
        cdef double x1, x2, x3, x4, x5, x6
        cdef ITYPE_t p, d1, d2
        cdef double A_p, q_p, A_d1, q_d1, A_d2, q_d2
        cdef double rho_
        cdef Py_ssize_t last_i, first_i
        cdef double w1, w2, x_i, W

        x1, x2, x3, x4, x5, x6 = x[0], x[1], x[2], x[3], x[4], x[5]
        p, d1, d2 = vessel_index_list[0], vessel_index_list[1], vessel_index_list[2]
        A_p, q_p, A_d1, q_d1, A_d2, q_d2 = u[0], u[1], u[2], u[3], u[4], u[5]
        rho_ = 1.

        last_i = self.vessels[p]._size - 1
        first_i = 0

        out[0] = -x2 + x4 + x6
        out[1] = (-0.5*rho_*(x2/x1)**2 - self.vessels[p].f_r0[last_i]*(1 -
                  sqrt(M_PI*self.vessels[p].r0[last_i]**2 / x1)) + 0.5*rho_*(x4/x3)**2 +
                  self.vessels[d1].f_r0[first_i]*(1 - sqrt(M_PI*self.vessels[d1].r0[0]**2 / x3)))
        out[2] = (-0.5*rho_*(x2/x1)**2 - self.vessels[p].f_r0[last_i]*(1 -
                  sqrt(M_PI*self.vessels[p].r0[last_i]**2 / x1)) + 0.5*rho_*(x6/x5)**2 +
                  self.vessels[d2].f_r0[0]*(1 - sqrt(M_PI*self.vessels[d2].r0[0]**2 / x5)))
        out[3] = (-x2/x1 - 4*self.pdes.compute_c_i(x1, index=last_i, vessel_index=p) +
                  (q_p/A_p + 4*self.pdes.compute_c_i(A_p, index=last_i-1, vessel_index=p)))
        out[4] = (-x4/x3 + 4*self.pdes.compute_c_i(x3, index=0, vessel_index=d1) +
                  (q_d1/A_d1 - 4*self.pdes.compute_c_i(A_d1, index=1, vessel_index=d1)))
        out[5] = (-x6/x5 + 4*self.pdes.compute_c_i(x5, index=0, vessel_index=d2) +
                  (q_d2/A_d2 - 4*self.pdes.compute_c_i(A_d2, index=1, vessel_index=d2)))

        return STATUS_OK

    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef ITYPE_t bifurcation_Jr(self, double[::1] x, ITYPE_t[::1] vessel_index_list, double[:, ::1] out):
        cdef double x1, x2, x3, x4, x5, x6, rho_
        cdef ITYPE_t p, d1, d2
        cdef Py_ssize_t last_i, first_i

        x1, x2, x3, x4, x5, x6 = x[0], x[1], x[2], x[3], x[4], x[5]
        p, d1, d2 = vessel_index_list[0], vessel_index_list[1], vessel_index_list[2]

        last_i = self.vessels[p]._size - 1
        first_i = 0

        rho_ = 1.
        out[0, 0] = 0.
        out[0, 1] = -1.
        out[0, 2] = 0.
        out[0, 3] = 1.
        out[0, 4] = 0.
        out[0, 5] = 1.
        out[1, 0] = rho_*(x2**2 / x1**3) - self.vessels[p].f_r0[last_i]*(M_PI*self.vessels[p].r0[last_i]**2)**0.5 / (2*x1**(3/2.))
        out[1, 1] = -rho_*x2 / x1**2
        out[1, 2] = -rho_*(x4**2 / x3**3) + self.vessels[d1].f_r0[0]*(M_PI*self.vessels[d1].r0[0]**2)**0.5 / (2*x3**(3/2.))
        out[1, 3] = rho_*x4 / x3**2
        out[1, 4] = 0.
        out[1, 5] = 0.
        out[2, 0] = rho_*(x2**2 / x1**3) - self.vessels[p].f_r0[last_i]*(M_PI*self.vessels[p].r0[last_i]**2)**0.5 / (2*x1**(3/2.))
        out[2, 1] = -rho_*x2 / x1**2
        out[2, 2] = 0.
        out[2, 3] = 0.
        out[2, 4] = -rho_*(x6**2 / x5**3) + self.vessels[d2].f_r0[0]*(M_PI*self.vessels[d2].r0[0]**2)**0.5 / (2*x5**(3/2.))
        out[2, 5] = rho_*x6 / x5**2
        out[3, 0] = (x2 / x1**2) - 4*(1/4.*sqrt((0.5/rho_)*self.vessels[p].f_r0[last_i]*sqrt(M_PI*self.vessels[p].r0[last_i]**2)))*x1**(-5/4.)
        out[3, 1] = -1./x1
        out[3, 2] = 0.
        out[3, 3] = 0.0
        out[3, 4] = 0.
        out[3, 5] = 0.
        out[4, 0] = 0.0
        out[4, 1] = 0.
        out[4, 2] = (x4 / x3**2) + 4*(1/4.*sqrt((0.5/rho_)*self.vessels[d1].f_r0[0]*sqrt(M_PI*self.vessels[d1].r0[0]**2)))*x3**(-5/4.)
        out[4, 3] = -1./x3
        out[4, 4] = 0.
        out[4, 5] = 0.
        out[5, 0] = 0.
        out[5, 1] = 0.
        out[5, 2] = 0.
        out[5, 3] = 0.
        out[5, 4] = (x6 / x5**2) + 4*(1/4.*sqrt((0.5/rho_)*self.vessels[d2].f_r0[0]*sqrt(M_PI*self.vessels[d2].r0[0]**2)))*x5**(-5/4.)
        out[5, 5] = -1./x5

        return STATUS_OK

    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cpdef ITYPE_t conjuction_R(self, double[::1] x, double[::1] u, double dt, ITYPE_t[::1] vessel_index_list, double[::1] out):
        cdef double x1, x2, x3, x4, rho_
        cdef double A_d1, q_d1, A_d2, q_d2
        cdef ITYPE_t d1, d2
        cdef Py_ssize_t last_i, first_i
        cdef double w1, w2, x_i, W

        x1, x2, x3, x4 = x[0], x[1], x[2], x[3]
        d1, d2 = vessel_index_list[0], vessel_index_list[1]

        #this is better to be improved and instead of shape - 1, last index should be shape
        last_i = self.vessels[d1]._size - 1
        first_i = 0

        A_d1, q_d1, A_d2, q_d2 = u[0], u[1], u[2], u[3]
        rho_ = 1.

        out[0] = -x2 + x4
        out[1] = (-0.5*rho_*(x2/x1)**2 - self.vessels[d1].f_r0[last_i]*(1 -
                  sqrt(M_PI*self.vessels[d1].r0[last_i]**2 / x1)) + 0.5*rho_*(x4/x3)**2 +
                  self.vessels[d2].f_r0[0]*(1 - sqrt(M_PI*self.vessels[d2].r0[0]**2 / x3)))
        out[2] = (-x2/x1 - 4*self.pdes.compute_c_i(x1, index=last_i, vessel_index=d1) +
                  (q_d1/A_d1 + 4*self.pdes.compute_c_i(A_d1, index=last_i-1, vessel_index=d1)))
        out[3] = (-x4/x3 + 4*self.pdes.compute_c_i(x3, index=0, vessel_index=d2) +
                  (q_d2/A_d2 - 4*self.pdes.compute_c_i(A_d2, index=1, vessel_index=d2)))

        return STATUS_OK

    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cpdef ITYPE_t conjuction_Jr(self, double[::1] x, ITYPE_t[::1] vessel_index_list, double[:,::1] out):
        cdef double x1, x2, x3, x4, rho_
        cdef ITYPE_t d1, d2
        cdef Py_ssize_t last_i, first_i

        x1, x2, x3, x4 = x[0], x[1], x[2], x[3]
        d1, d2 = vessel_index_list[0], vessel_index_list[1]
        rho_ = 1.
        last_i = self.vessels[d1]._size - 1

        out[0, 0] = 0.
        out[0, 1] = -1.
        out[0, 2] = 0.
        out[0, 3] = 1.
        out[1, 0] = rho_*(x2**2 / x1**3) - self.vessels[d1].f_r0[last_i]*(M_PI*self.vessels[d1].r0[last_i]**2)**0.5 / (2*x1**(3/2.))
        out[1, 1] = -rho_*x2 / x1**2
        out[1, 2] = -rho_*(x4**2 / x3**3) + self.vessels[d2].f_r0[0]*(M_PI*self.vessels[d2].r0[0]**2)**0.5 / (2*x3**(3/2.))
        out[1, 3] = rho_*x4 / x3**2
        out[2, 0] = (x2 / x1**2) - 4*(1/4.*sqrt((0.5/rho_)*self.vessels[d1].f_r0[last_i]*sqrt(M_PI*self.vessels[d1].r0[last_i]**2)))*x1**(-5/4.)
        out[2, 1] = -1./x1
        out[2, 2] = 0.
        out[2, 3] = 0.
        out[3, 0] = 0.
        out[3, 1] = 0.
        out[3, 2] = (x4 / x3**2) + 4*(1/4.*sqrt((0.5/rho_)*self.vessels[d2].f_r0[0]*sqrt(M_PI*self.vessels[d2].r0[0]**2)))*x3**(-5/4.)
        out[3, 3] = -1./x3

        return STATUS_OK


cdef class BCsD(cBCs):
    cdef void U_0(self, double** u, Py_ssize_t clms, double t, double dx, double dt, ITYPE_t vessel_index, double** out):
        theta = dt / dx
        dt2 = dt / 2.
        p_pres = self.inlet_fun(t)

        A0 = np.pi*self.mesh.vessels[vessel_index].r0[0]**2.0
        f = self.mesh.vessels[vessel_index].f_r0[0]

        A = A0 / (1 - (p_pres /f))**2.0
        # A = ((p_pres/(f*np.sqrt(A0))) + np.sqrt(A0))**2
        W2_0 = u[1][0]/u[0][0] - 4*self.pdes.compute_c_i(A, index=0, vessel_index=vessel_index)
        W2_1 = u[1][1]/u[0][1] - 4*self.pdes.compute_c_i(u[0][1], index=1, vessel_index=vessel_index)
        l2 = u[1][0]/u[0][0] - self.pdes.compute_c_i(A, index=0, vessel_index=vessel_index)
        x = fabs(-l2*dt)
        W2 = linear_extrapolation(x, 0.0, dx, W2_0, W2_1)
        q = A*(W2 + 4*self.pdes.compute_c_i(A, index=0, vessel_index=vessel_index))

        out[0][0] = A
        out[1][0] = q


cdef class BCsADAN56(cBCsWat):

    cdef void U_0(self, double** u, Py_ssize_t clms, double t, double dx, double dt, ITYPE_t vessel_index, double** out):
        cdef double theta, dt2, q_1, A, q_pres

        theta = dt/dx
        dt2 = 0.5*dt

        q_pres = self.inlet_fun(t)

        q_1 = out[1][1]

        # A[n+1, 0]
        A = u[0][0] - (2.*theta)*(q_1 - q_pres)

        out[0][0] = A
        out[1][0] = q_pres

    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef ITYPE_t U_L(self, double** u, Py_ssize_t clms, double t, double dx, double dt,
                     ITYPE_t vessel_index, double** out) nogil except *:
        cdef double theta, dt2, out_a, out_q, p0, p_out, x
        cdef Py_ssize_t rws, i, j
        cdef double q_m1
        rws = 2
        theta = dt/dx
        dt2 = dt/2.
        q_m1 = out[1][clms-2]

        p_out = self.pdes.pressure_i(u[0][clms-1], clms-1, vessel_index)
        p0 = self.pdes.pressure_i(u[0][clms-1], clms-1, vessel_index)

        cdef Py_ssize_t k = 0
        cdef double R1, Rt, Ct, p_old

        R1 = self.vessels[vessel_index].RLC[0]
        Rt = self.vessels[vessel_index].RLC[1]
        Ct = self.vessels[vessel_index].RLC[2]

        x = (dt / (R1*Rt*Ct))
        while k < 1000:
            p_old = p0
            out_q = (x*p_out - x*(R1 + Rt)*u[1][clms-1] +
                     (p0 - p_out)/R1 + u[1][clms-1])

            out_a = u[0][clms-1] - theta*(out_q - q_m1)
            p0 = self.pdes.pressure_i(out_a, index=clms-1, vessel_index=vessel_index)
            if fabs(p_old - p0) < 1e-7:
                break
            k += 1
        out[0][clms-1] = out_a
        out[1][clms-1] = out_q

        return STATUS_OK

cdef class cBCsHandModelNonReflBcs(cBCsWat):
    cdef void U_0(self, double** u, Py_ssize_t clms, double t, double dx, double dt, ITYPE_t vessel_index, double** out):
        theta = dt / dx
        dt2 = dt / 2.

        p_pres = self.inlet_fun(t)

        A0 = np.pi*self.mesh.vessels[vessel_index].r0[0]**2.0
        A1 = np.pi*self.mesh.vessels[vessel_index].r0[1]**2.0
        f = self.mesh.vessels[vessel_index].f_r0[0]

        A = A0*((p_pres /f) + 1)**2.0

        W2_const = -4.0*self.pdes.compute_c_i(A1, index=1, vessel_index=vessel_index)

        W2 = u[1][1]/u[0][1] - 4*self.pdes.compute_c_i(u[0][1], index=1, vessel_index=vessel_index)

        W1 = W2_const + 8*np.sqrt(f/(2*self._rho))*(((p_pres/f) + 1)**(1/2.))

        A_out = A0*(((W1 - W2)/(4.))**4)*(((self._rho)/(2*f))**2)
        q_out = 0.5*A_out*(W1 + W2)

        out[0][0] = A_out
        out[1][0] = q_out


cdef class cFDMSolver:
    """
    Base class for finite-difference computing schemes.
    """
    cdef Pool mem
    cdef object mesh
    cdef cPDEm _pdes
    cdef cBCs _bcs
    cdef Py_ssize_t _no_cycles
    cdef double _dt
    cdef Py_ssize_t _Nx
    cdef Py_ssize_t _Nt
    cdef double[::1] _t
    cdef double _T
    cdef ITYPE_t[::1] U0_arr
    cdef ITYPE_t[::1] UL_arr
    cdef ITYPE_t[::1] _It
    cdef ITYPE_t[:, ::1] UBif_arr
    cdef ITYPE_t[:, ::1] UConj_arr
    cdef double _C

    def __init__(self, cBCs bcs):
        self.mem = Pool()
        self._bcs = bcs
        self._pdes = self._bcs.get_pdes()
        self.mesh = self._bcs.get_pdes().mesh

    @property
    def mesh(self):
        return self.mesh

    cpdef void set_T(self, double dt, double T, Py_ssize_t no_cycles):
        self._no_cycles = no_cycles
        self._T = T
        self._dt = dt
        self._Nt = int(round(self._T / self._dt))
        self._t = np.linspace(0.0, self._Nt*self._dt, self._Nt + 1)
        self._pdes.set_boundary_layer_th(T, no_cycles)

    cpdef void set_BC(self, ITYPE_t[::1] U0_array,
                      ITYPE_t[::1] UL_array, ITYPE_t[:, ::1] UBif_array,
                      ITYPE_t[:,::1] UConj_array):
        self.U0_arr = U0_array
        self.UL_arr = UL_array
        self.UBif_arr = UBif_array
        self.UConj_arr = UConj_array

    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef CFL cfl_i(self, double[::1] uin, double dx, double dt, ITYPE_t index, ITYPE_t vessel_index):
        cdef double _A, _q, theta, c, u, ls, rs, comp_elem
        cdef double[::1] W = np.zeros(2, dtype=np.float)
        cdef ITYPE_t i
        _A = uin[0]
        _q = uin[1]
        theta = dt / dx
        c = self._pdes.compute_c_i(_A, index, vessel_index)
        u = _q / _A
        W[0] = abs(u + c)
        W[1] = abs(u - c)
        if ((1./W[0]) <= (1./W[1])) == True:
            comp_elem = 1./W[0]
        else:
            comp_elem = 1./W[1]
        if (theta < comp_elem) == True:
            return CFL.CFL_STATUS_OK
        else:
            return CFL.CFL_STATUS_ERROR

    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef CFL cfl_condition(self, double** uin, Py_ssize_t clms, double dx, double dt, ITYPE_t vessel_index):
        cdef double _A, _q, theta, c, u, ls, rs, comp_elem
        cdef double[::1] W = np.zeros(2, dtype=np.float)
        cdef Py_ssize_t i
        cdef Py_ssize_t size
        cdef double min_value

        size = clms
        W = np.zeros(size, dtype=np.float)
        theta = dt / dx

        c = self._pdes.compute_c_i(uin[0][0], 0, vessel_index)
        W[0] = fabs(1.0 / ((uin[1][0]/uin[0][0]) + c))
        min_value = W[0]

        for i in range(1, size):
            c = self._pdes.compute_c_i(uin[0][i], i, vessel_index)
            W[i] = fabs(1.0 / ((uin[1][i]/uin[0][i]) + c))
            if W[i] < min_value:
                min_value = W[i]

        if (theta < min_value) == True:
            return CFL.CFL_STATUS_OK
        else:
            return CFL.CFL_STATUS_ERROR

    cdef void initialise_solution_vector(self, vector[U] &vectorU, ITYPE_t length):
        raise NotImplementedError("This method has not been implemented!")


    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.initializedcheck(False)
    @cython.cdivision(True)
    cpdef SOLVER solve(self, user_action=None,
              version="vectorised"):

        cdef Pool mem = Pool()
        cdef Py_ssize_t length, len_x, i, j, n, i_bcs
        cdef double[::1] getarr = np.zeros(2, np.float)
        # array to store the wave speeds
        cdef double[::1] c_i
        cdef CFL res_cfl
        cdef SOLVER res_solver
        
        length = len(self.mesh.vessels)

        cdef Py_ssize_t siz_i, i_save
        cdef vector[U] v_U
        cdef vector[vector[Solution]] vSol
        v_U.resize(length)
        

        self.initialise_solution_vector(v_U, length)

        import time
        t0 = time.perf_counter()  # CPU time measurement

        # --- Valid indices for space and time mesh ---
        self._It = np.arange(0, self._Nt + 1, dtype=np.int64)

        self._dt = self._t[1] - self._t[0]

        print("Solver set to dt=%0.9f" % self._dt)

        if user_action is not None:
            i_save = 0
            vSol.push_back(vector[Solution]())
            vSol[i_save].resize(length)

        # check with nogil here
        for i in range(length):
            len_x = self.mesh.vessels[i].x.shape[0]
            if user_action is not None:
                vSol[i_save][i].n = 0
                vSol[i_save][i].size = 3*len_x
                vSol[i_save][i].u_store = <double*>mem.alloc(3*len_x, sizeof(double)) 
            for j in range(len_x):
                self._bcs.I(self.mesh.vessels[i].x[j], j, i, getarr)
                v_U[i].u_n[0][j] = getarr[0]
                v_U[i].u_n[1][j] = getarr[1]

                # check cfl condition
                res_cfl = self.cfl_i(getarr, v_U[i].dx, self._dt, j, i)
                if res_cfl == CFL.CFL_STATUS_ERROR:
                    return SOLVER.SOLVER_STATUS_ERROR

                if user_action is not None:
                    p = self._pdes.pressure_i(getarr[0], j, vessel_index=i)
                    
                    vSol[i_save][i].u_store[j] = v_U[i].u_n[0][j]
                    vSol[i_save][i].u_store[j + len_x] = v_U[i].u_n[1][j]
                    vSol[i_save][i].u_store[j + 2*len_x] = p


#         # TIME LOOP
        cdef ITYPE_t index, siz_U0, k, siz_UL, siz_UBif, siz_UConj, par, d1, d2
        cdef double[::1] u_nmp1, x_n, u_nmp1_conj, x_n_conj
        u_nmp1 = np.zeros(6, np.float)
        x_n = np.zeros(6, np.float)
        u_nmp1_conj = np.zeros(4, np.float)
        x_n_conj = np.zeros(4, np.float)
        cdef ITYPE_t last_index, first_index, iter_count
        cdef ITYPE_t[::1] bif_list = np.zeros(3, dtype=np.int64)
        cdef ITYPE_t[::1] conj_list = np.zeros(2, dtype=np.int64)
        siz_U0 = self.U0_arr.shape[0]
        siz_UL = self.UL_arr.shape[0]
        siz_UBif = self.UBif_arr.shape[0]
        siz_UConj = self.UConj_arr.shape[0]

        #range [1 , self._Nt + 1]
        for n in range(1, self._Nt + 1):
            # this needs to be parallel
            # with nogil, parallel(num_threads=8):
            for i in prange(length, nogil=True):
            # for i in range(length):
                self.advance(v_U[i], self._dt, v_U[i].dx, i)

#             # Insert boundary conditions
            for index in range(siz_U0):
                k = self.U0_arr[index]
                self._bcs.U_0(v_U[k].u_n, v_U[k].x_size, self._t[n], v_U[k].dx, self._dt, k, v_U[k].u)

            # with nogil, parallel(num_threads=8):
            for index in prange(siz_UL, nogil=True):
                self._bcs.U_L(v_U[self.UL_arr[index]].u_n, v_U[self.UL_arr[index]].x_size, self._t[n], v_U[self.UL_arr[index]].dx,
                              self._dt, self.UL_arr[index], v_U[self.UL_arr[index]].u)

            #parallel siz_UBif
            for index in range(siz_UBif):
            # for index in prange(siz_UBif, nogil=True):
                par = self.UBif_arr[index, 0]
                d1 = self.UBif_arr[index, 1]
                d2 = self.UBif_arr[index, 2]
                bif_list[0] = par
                bif_list[1] = d1
                bif_list[2] = d2
                u_nmp1[0], u_nmp1[1], u_nmp1[2], u_nmp1[3], u_nmp1[4], u_nmp1[5] = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
                x_n[0], x_n[1], x_n[2], x_n[3], x_n[4], x_n[5] = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
                last_index = v_U[par].x_size
                first_index = 0
                u_nmp1[0] = v_U[par].u[0][last_index - 2]
                u_nmp1[1] = v_U[par].u[1][last_index - 2]
                # print v_U[par].u[1, last_index - 2]
                x_n[0] = v_U[par].u_n[0][last_index - 1]
                x_n[1] = v_U[par].u_n[1][last_index - 1]
                u_nmp1[2] = v_U[d1].u[0][first_index + 1]
                u_nmp1[3] = v_U[d1].u[1][first_index + 1]
                x_n[2] = v_U[d1].u_n[0][first_index]
                x_n[3] = v_U[d1].u_n[1][first_index]
                u_nmp1[4] = v_U[d2].u[0][first_index + 1]
                u_nmp1[5] = v_U[d2].u[1][first_index + 1]
                x_n[4] = v_U[d2].u_n[0][first_index]
                x_n[5] = v_U[d2].u_n[1][first_index]


                # return SOLVER.SOLVER_STATUS_ERROR
                iter_count = Newton_Raphson_Sys_bif(self._bcs, x_n, u_nmp1, self._dt,
                                                    bif_list, 1.0e-9, 100)

                # print iter_count
                v_U[par].u[0][last_index - 1] = x_n[0]
                v_U[par].u[1][last_index - 1] = x_n[1]
                v_U[d1].u[0][0] = x_n[2]
                v_U[d1].u[1][0] = x_n[3]
                v_U[d2].u[0][0] = x_n[4]
                v_U[d2].u[1][0] = x_n[5]

            for index in range(siz_UConj):
                d1 = self.UConj_arr[index, 0]
                d2 = self.UConj_arr[index, 1]
                conj_list[0] = d1
                conj_list[1] = d2
                u_nmp1_conj[0], u_nmp1_conj[1], u_nmp1_conj[2], u_nmp1_conj[3] = 0.0, 0.0, 0.0, 0.0
                x_n_conj[0], x_n_conj[1], x_n_conj[2], x_n_conj[3] = 0.0, 0.0, 0.0, 0.0
                last_index = v_U[d1].x_size
                first_index = 0
                u_nmp1_conj[0] = v_U[d1].u[0][last_index - 2]
                u_nmp1_conj[1] = v_U[d1].u[1][last_index - 2]
                x_n_conj[0] = v_U[d1].u_n[0][last_index - 1]
                x_n_conj[1] = v_U[d1].u_n[1][last_index - 1]
                u_nmp1_conj[2] = v_U[d2].u[0][first_index + 1]
                u_nmp1_conj[3] = v_U[d2].u[1][first_index + 1]
                x_n_conj[2] = v_U[d2].u_n[0][first_index]
                x_n_conj[3] = v_U[d2].u_n[1][first_index]

                iter_count = Newton_Raphson_Sys_conj(self._bcs, x_n_conj, u_nmp1_conj, self._dt,
                                                     conj_list, 1.0e-9, 100)


                v_U[d1].u[0][last_index - 1] = x_n_conj[0]
                v_U[d1].u[1][last_index - 1] = x_n_conj[1]
                v_U[d2].u[0][0] = x_n_conj[2]
                v_U[d2].u[1][0] = x_n_conj[3]

            if (user_action is not None) and (n % user_action.skip_frame == 0):
                i_save += 1
                vSol.push_back(vector[Solution]())
                vSol[i_save].resize(length)

            for i in range(length):
                if (user_action is not None) and (n % user_action.skip_frame == 0):
                    # # check cfl condition
                    res_cfl = self.cfl_condition(v_U[i].u, v_U[i].x_size, v_U[i].dx, self._dt, i)
                    if res_cfl == CFL.CFL_STATUS_ERROR:
                        print("Solver failed in vessel %d in time increment %d  spatial %d and lenght=%d" % (i, n, j, v_U[i].x_size))
                        return SOLVER.SOLVER_STATUS_ERROR
                    else:
                        vSol[i_save][i].n = n                        
                        vSol[i_save][i].size = 3*v_U[i].x_size
                        vSol[i_save][i].u_store = <double*>mem.alloc(vSol[i_save][i].size, sizeof(double)) 
                        for j in range(v_U[i].x_size):
                            getarr[0] = v_U[i].u[0][j]
                            getarr[1] = v_U[i].u[1][j]

                            p = self._pdes.pressure_i(getarr[0], j, vessel_index=i)

                            vSol[i_save][i].u_store[j] = v_U[i].u[0][j]
                            vSol[i_save][i].u_store[j + v_U[i].x_size] = v_U[i].u[1][j]
                            vSol[i_save][i].u_store[j + 2*v_U[i].x_size] = p

                for j in range(v_U[i].x_size):
                    v_U[i].u_n[0][j] = v_U[i].u[0][j]
                    v_U[i].u_n[1][j] = v_U[i].u[1][j]


            printf("\r%6.2f%%", (100.0 * self._t[n] / self._T))

        cdef Py_ssize_t size_solution = vSol.size()
        cdef Py_ssize_t id, jd
        cdef double[:, ::1] u_store
        if user_action is not None:

            print("Size of solution is:", size_solution)
            print("Writing solution to file...")
            
            for id in range(size_solution):
                for jd in range(length):
                    u_store = <double[:3, :vSol[id][jd].size/3]>(vSol[id][jd].u_store)
                    user_action(u_store, self.mesh.vessels[jd].x, self._t, vSol[id][jd].n, PRINT_STATUS, WRITE_STATUS, jd)
            print("Writing solution to file completed!")

        return SOLVER.SOLVER_STATUS_OK


    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void advance(self, U vec, double dt, double dx, ITYPE_t vessel_index) nogil:
        cdef Py_ssize_t siz
        cdef Py_ssize_t i, j

        cdef double dt2, theta

        siz = vec.x_size

        dt2 = dt/2.

        theta = dt / dx

        self._pdes.flux(vec.u_n, vec.x_size, 0, vessel_index, vec.F_)
        self._pdes.source(vec.u_n, vec.x_size, 0, vessel_index, vec.S_)

#         # U[n+1/2, i-1/2]
#         with nogil, parallel(num_threads=4):
        for i in range(1, vec.x_size-1):
            vec.u_nh_mh[0][i] = ((vec.u_n[0][i] + vec.u_n[0][i-1]) / 2. -
                                 0.5 * (theta) * (vec.F_[0][i] - vec.F_[0][i-1]) +
                                 0.5 * dt2 * (vec.S_[0][i] + vec.S_[0][i-1]))
            vec.u_nh_mh[1][i] = ((vec.u_n[1][i] + vec.u_n[1][i-1]) / 2. -
                                 0.5 * (theta) * (vec.F_[1][i] - vec.F_[1][i-1]) +
                                 0.5 * dt2 * (vec.S_[1][i] + vec.S_[1][i-1]))

            vec.u_nh_ph[0][i] = ((vec.u_n[0][i+1] + vec.u_n[0][i]) / 2. -
                                 0.5 * (theta) * (vec.F_[0][i+1] - vec.F_[0][i]) +
                                 0.5 * dt2 * (vec.S_[0][i+1] + vec.S_[0][i]))
            vec.u_nh_ph[1][i] = ((vec.u_n[1][i+1] + vec.u_n[1][i]) / 2. -
                                 0.5 * (theta) * (vec.F_[1][i+1] - vec.F_[1][i]) +
                                 0.5 * dt2 * (vec.S_[1][i+1] + vec.S_[1][i]))


        # this is a test, because we need to fill the edges with ghost
        vec.u_nh_mh[0][0] = vec.u_n[0][0]
        vec.u_nh_mh[1][0] = vec.u_n[1][0]
        vec.u_nh_mh[0][siz-1] = vec.u_n[0][siz-1]
        vec.u_nh_mh[1][siz-1] = vec.u_n[1][siz-1]
        vec.u_nh_ph[0][0] = vec.u_n[0][0]
        vec.u_nh_ph[1][0] = vec.u_n[1][0]
        vec.u_nh_ph[0][siz-1] = vec.u_n[0][siz-1]
        vec.u_nh_ph[1][siz-1] = vec.u_n[1][siz-1]


        self._pdes.flux(vec.u_nh_mh, vec.x_size, -1, vessel_index, vec.F_nh_mh)

        self._pdes.source(vec.u_nh_mh, vec.x_size, -1, vessel_index, vec.S_nh_mh)

        self._pdes.flux(vec.u_nh_ph, vec.x_size, 1, vessel_index, vec.F_nh_ph)

        self._pdes.source(vec.u_nh_ph, vec.x_size, 1, vessel_index, vec.S_nh_ph)


        for i in range(1, siz-1):
            vec.u[0][i] = (vec.u_n[0][i] - theta * (vec.F_nh_ph[0][i] - vec.F_nh_mh[0][i]) +
                           dt2*(vec.S_nh_ph[0][i] + vec.S_nh_mh[0][i]))
            vec.u[1][i] = (vec.u_n[1][i] - theta * (vec.F_nh_ph[1][i] - vec.F_nh_mh[1][i]) +
                           dt2*(vec.S_nh_ph[1][i] + vec.S_nh_mh[1][i]))


cdef class cLaxWendroffSolver(cFDMSolver):
    """

    Class with 2 step Lax-Wendroff scheme.

    .. math::

        U_i^{n+1} = U_i^{n} - \\frac{\\Delta t}{\\Delta x} \\left( F_{i + \\frac{1}{2}}^{n + \\frac{1}{2}}
         - F_{i - \\frac{1}{2}}^{n + \\frac{1}{2}} \\right) +
          \\frac{\\Delta t}{2}\\left( S_{i + \\frac{1}{2}}^{n + \\frac{1}{2}}
           + S_{i - \\frac{1}{2}}^{n + \\frac{1}{2}} \\right)

    where the intermediate values calculated as

    .. math::

        U_j^{n+\\frac{1}{2}} = \\frac{U_{j+1/2}^n + U_{j-1/2}^n}{2}
         - \\frac{\\Delta t}{2 \\Delta x}\\left( F_{j+1/2} - F_{j-1/2} \\right)
          + \\frac{\\Delta t}{4}\\left( S_{j+1/2} + S_{j-1/2}  \\right)

    """
    cdef void initialise_solution_vector(self, vector[U] &vectorU, ITYPE_t length):
        cdef Py_ssize_t i, j, k
        for i in range(length):
            siz_i = self.mesh.vessels[i].x.shape[0]
            vectorU[i].dx = self.mesh.vessels[i].dx
            vectorU[i].x_size = siz_i
            vectorU[i].u = <double**>self.mem.alloc(2, sizeof(double*))
            vectorU[i].u_n =  <double**>self.mem.alloc(2, sizeof(double*))
            vectorU[i].F_ =  <double**>self.mem.alloc(2, sizeof(double*))
            vectorU[i].S_ =  <double**>self.mem.alloc(2, sizeof(double*))
            vectorU[i].u_nh_mh =  <double**>self.mem.alloc(2, sizeof(double*))
            vectorU[i].u_nh_ph =  <double**>self.mem.alloc(2, sizeof(double*))
            vectorU[i].F_nh_mh =  <double**>self.mem.alloc(2, sizeof(double*))
            vectorU[i].S_nh_mh =  <double**>self.mem.alloc(2, sizeof(double*))
            vectorU[i].F_nh_ph =  <double**>self.mem.alloc(2, sizeof(double*))
            vectorU[i].S_nh_ph =  <double**>self.mem.alloc(2, sizeof(double*))
            for j in range(2):
                vectorU[i].u[j] = <double*>self.mem.alloc(siz_i, sizeof(double))
                vectorU[i].u_n[j] =  <double*>self.mem.alloc(siz_i, sizeof(double))
                vectorU[i].F_[j] =  <double*>self.mem.alloc(siz_i, sizeof(double))
                vectorU[i].S_[j] =  <double*>self.mem.alloc(siz_i, sizeof(double))
                vectorU[i].u_nh_mh[j] =  <double*>self.mem.alloc(siz_i, sizeof(double))
                vectorU[i].u_nh_ph[j] =  <double*>self.mem.alloc(siz_i, sizeof(double))
                vectorU[i].F_nh_mh[j] =  <double*>self.mem.alloc(siz_i, sizeof(double))
                vectorU[i].S_nh_mh[j] =  <double*>self.mem.alloc(siz_i, sizeof(double))
                vectorU[i].F_nh_ph[j] =  <double*>self.mem.alloc(siz_i, sizeof(double))
                vectorU[i].S_nh_ph[j] =  <double*>self.mem.alloc(siz_i, sizeof(double))


    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void advance(self, U vec, double dt, double dx, ITYPE_t vessel_index) nogil:
        cdef Py_ssize_t siz
        cdef Py_ssize_t i, j

        cdef double dt2, theta

        siz = vec.x_size

        dt2 = dt/2.

        theta = dt / dx

        self._pdes.flux(vec.u_n, siz, 0, vessel_index, vec.F_)
        self._pdes.source(vec.u_n, siz, 0, vessel_index, vec.S_)

#         # U[n+1/2, i-1/2]
#         with nogil, parallel(num_threads=4):
        for i in range(1, vec.x_size-1):
            vec.u_nh_mh[0][i] = ((vec.u_n[0][i] + vec.u_n[0][i-1]) / 2. -
                                 0.5 * (theta) * (vec.F_[0][i] - vec.F_[0][i-1]) +
                                 0.5 * dt2 * (vec.S_[0][i] + vec.S_[0][i-1]))
            vec.u_nh_mh[1][i] = ((vec.u_n[1][i] + vec.u_n[1][i-1]) / 2. -
                                 0.5 * (theta) * (vec.F_[1][i] - vec.F_[1][i-1]) +
                                 0.5 * dt2 * (vec.S_[1][i] + vec.S_[1][i-1]))

            vec.u_nh_ph[0][i] = ((vec.u_n[0][i+1] + vec.u_n[0][i]) / 2. -
                                 0.5 * (theta) * (vec.F_[0][i+1] - vec.F_[0][i]) +
                                 0.5 * dt2 * (vec.S_[0][i+1] + vec.S_[0][i]))
            vec.u_nh_ph[1][i] = ((vec.u_n[1][i+1] + vec.u_n[1][i]) / 2. -
                                 0.5 * (theta) * (vec.F_[1][i+1] - vec.F_[1][i]) +
                                 0.5 * dt2 * (vec.S_[1][i+1] + vec.S_[1][i]))


        # this is a test, because we need to fill the edges with ghost
        vec.u_nh_mh[0][0] = vec.u_n[0][0]
        vec.u_nh_mh[1][0] = vec.u_n[1][0]
        vec.u_nh_mh[0][siz-1] = vec.u_n[0][siz-1]
        vec.u_nh_mh[1][siz-1] = vec.u_n[1][siz-1]
        vec.u_nh_ph[0][0] = vec.u_n[0][0]
        vec.u_nh_ph[1][0] = vec.u_n[1][0]
        vec.u_nh_ph[0][siz-1] = vec.u_n[0][siz-1]
        vec.u_nh_ph[1][siz-1] = vec.u_n[1][siz-1]


        self._pdes.flux(vec.u_nh_mh, siz, -1, vessel_index, vec.F_nh_mh)

        self._pdes.source(vec.u_nh_mh, siz, -1, vessel_index, vec.S_nh_mh)

        self._pdes.flux(vec.u_nh_ph, siz, 1, vessel_index, vec.F_nh_ph)

        self._pdes.source(vec.u_nh_ph, siz, 1, vessel_index, vec.S_nh_ph)

        for i in range(1, siz-1):
            vec.u[0][i] = (vec.u_n[0][i] - theta * (vec.F_nh_ph[0][i] - vec.F_nh_mh[0][i]) +
                           dt2*(vec.S_nh_ph[0][i] + vec.S_nh_mh[0][i]))
            vec.u[1][i] = (vec.u_n[1][i] - theta * (vec.F_nh_ph[1][i] - vec.F_nh_mh[1][i]) +
                           dt2*(vec.S_nh_ph[1][i] + vec.S_nh_mh[1][i]))

cdef class cMacCormackSolver(cFDMSolver):
    """
    Class with the MacCormack scheme implementation.

    .. math::

        U_i^{\\star} = U_i^n - \\frac{\\Delta t}{\\Delta x}\\left( F_{i+1}^n -
         F_i^n \\right) + \\Delta t S_i^n



    and

    .. math::

        U_i^{n+1} = \\frac{1}{2}\\left( U_i^n + U_i^{\\star} \\right)
         - \\frac{\\Delta t}{2 \\Delta x}\\left( F_i^{\\star}
          - F_{i-1}^{\\star} \\right) + \\frac{\\Delta t}{2}S_i^{\\star}

    """
    cdef void initialise_solution_vector(self, vector[U] &vectorU, ITYPE_t length):
        cdef Py_ssize_t i, j
        for i in range(length):
            siz_i = self.mesh.vessels[i].x.shape[0]
            vectorU[i].x_size = siz_i
            vectorU[i].dx = self.mesh.vessels[i].dx
            vectorU[i].u = <double**>self.mem.alloc(2, sizeof(double*))
            vectorU[i].u_n = <double**>self.mem.alloc(2, sizeof(double*))
            vectorU[i].u_star = <double**>self.mem.alloc(2, sizeof(double*))
            vectorU[i].F_ = <double**>self.mem.alloc(2, sizeof(double*))
            vectorU[i].S_ = <double**>self.mem.alloc(2, sizeof(double*))
            vectorU[i].F_star = <double**>self.mem.alloc(2, sizeof(double*))
            vectorU[i].S_star = <double**>self.mem.alloc(2, sizeof(double*))
            for j in range(2):
                vectorU[i].u[j] = <double*>self.mem.alloc(siz_i, sizeof(double))
                vectorU[i].u_n[j] = <double*>self.mem.alloc(siz_i, sizeof(double))
                vectorU[i].u_star[j] = <double*>self.mem.alloc(siz_i, sizeof(double))
                vectorU[i].F_[j] = <double*>self.mem.alloc(siz_i, sizeof(double))
                vectorU[i].S_[j] = <double*>self.mem.alloc(siz_i, sizeof(double))
                vectorU[i].F_star[j] = <double*>self.mem.alloc(siz_i, sizeof(double))
                vectorU[i].S_star[j] = <double*>self.mem.alloc(siz_i, sizeof(double))



    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void advance(self, U vec, double dt, double dx, ITYPE_t vessel_index) nogil:
        cdef Py_ssize_t siz
        cdef Py_ssize_t i, j

        cdef double dt2, theta

        siz = vec.x_size

        dt2 = 0.5*dt

        theta = dt / dx

        self._pdes.flux(vec.u_n, siz, 0, vessel_index, vec.F_)
        self._pdes.source(vec.u_n, siz, 0, vessel_index, vec.S_)

#         # U[n+1/2, i-1/2]
#         with nogil, parallel(num_threads=4):
        for i in range(0, siz-1):
            vec.u_star[0][i] = (vec.u_n[0][i] - theta*(vec.F_[0][i+1] - vec.F_[0][i]) +
                                dt*vec.S_[0][i])
            vec.u_star[1][i] = (vec.u_n[1][i] - theta*(vec.F_[1][i+1] - vec.F_[1][i]) +
                                dt*vec.S_[1][i])

        vec.u_star[0][siz-1] = vec.u_n[0][siz-1]
        vec.u_star[1][siz-1] = vec.u_n[1][siz-1]

        vec.u_star[0][siz-1] = (vec.u_n[0][siz-1] - theta*(vec.F_[0][siz-2] - vec.F_[0][siz-1]) +
                                dt*vec.S_[0][siz-1])
        vec.u_star[1][siz-1] = (vec.u_n[1][siz-1] - theta*(vec.F_[1][siz-2] - vec.F_[1][siz-1]) +
                                dt*vec.S_[1][siz-1])


        self._pdes.flux(vec.u_star, siz, 0, vessel_index, vec.F_star)
        self._pdes.source(vec.u_star, siz, 0, vessel_index, vec.S_star)

        for i in range(1, siz-1):
            vec.u[0][i] = (0.5*(vec.u_n[0][i] + vec.u_star[0][i]) - 0.5*theta*(vec.F_star[0][i] -
                           vec.F_star[0][i-1]) + 0.5*dt * vec.S_star[0][i])
            vec.u[1][i] = (0.5*(vec.u_n[1][i] + vec.u_star[1][i]) - 0.5*theta*(vec.F_star[1][i] -
                           vec.F_star[1][i-1]) + 0.5*dt * vec.S_star[1][i])


cdef class cMacCormackGodunovSplitSolver(cMacCormackSolver):
    """
    Solve U_t + F_x = S
    """
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.initializedcheck(False)
    @cython.cdivision(True)
    cdef ITYPE_t _advance_Godunov_splitting(self, cParabolicSolObj parStruct, double** u, Py_ssize_t N_n, ITYPE_t n, double dt,
                                             double dx, ITYPE_t vessel_index, double theta):
        cdef Pool mem = Pool()
        cdef Py_ssize_t i
        cdef ITYPE_t res
        cdef double F, gamma
        cdef double* C_v = <double*>mem.alloc(N_n, sizeof(double))
        cdef double* a = <double*>mem.alloc(N_n, sizeof(double))
        cdef double* q
        cdef double* d
        
        F = dt / (dx * dx)

        for i in range(N_n):
            a[i] = u[0][i]

        res = self._pdes.CV_f(a, N_n, vessel_index, C_v)

        for i in range(N_n-2):
            parStruct.a_ph[i] = 0.5 * (C_v[i+1] + C_v[i+2])
        # print(a_ph)
            parStruct.a_mh[i] = 0.5 * (C_v[i+1] + C_v[i])
            parStruct.b[i] = (1 + (u[0][i+1]/self._pdes._rho)*theta * F *
                              (parStruct.a_ph[i] + parStruct.a_mh[i]))

        for i in range(N_n-3):
            parStruct.a[i] = -(u[0][i+2]/self._pdes._rho)*theta * F * parStruct.a_mh[i+1]
            parStruct.c[i] = -(u[0][i+1]/self._pdes._rho)*theta * F * parStruct.a_ph[i]

        # # ---- BCd Neuman ----- ##
        parStruct.b[0] = (1 + (u[0][1]/self._pdes._rho)*theta * F *
                          (parStruct.a_ph[0] + parStruct.a_mh[0])) -\
                         (u[0][1]/self._pdes._rho)*theta * F * parStruct.a_mh[0]
        # here the value for m, m-1
        parStruct.b[N_n-3] = (1 + (u[0][N_n-2]/self._pdes._rho)*theta * F *
                              (parStruct.a_ph[N_n-3] + parStruct.a_mh[N_n-3])) -\
                             (u[0][N_n-2]/self._pdes._rho)*theta * F * parStruct.a_ph[N_n-3]

        gamma = (1 - theta)

        for i in range(1, N_n-1):
            parStruct.d[i] = ((u[0][i]/self._pdes._rho)*gamma * F * parStruct.a_mh[i-1] * u[1][i-1] +
                              (1 - (u[0][i]/self._pdes._rho)*gamma * F * (parStruct.a_mh[i-1] +
                              parStruct.a_ph[i-1])) * u[1][i] +
                              (u[0][i]/self._pdes._rho)*gamma * F * parStruct.a_ph[i-1] * u[1][i+1])
        

        q = <double*>mem.alloc(N_n-2, sizeof(double))
        d = <double*>mem.alloc(N_n-2, sizeof(double))
        for i in range(1, N_n - 1):
            q[i-1] = u[1][i]
            d[i-1] = parStruct.d[i]

        res = tdma(parStruct.a, parStruct.b, parStruct.c, d, q, N_n-2)
        for i in range(1, N_n - 1):
            u[1][i] = q[i-1]


    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.initializedcheck(False)
    @cython.cdivision(True)
    cpdef SOLVER solve(self, user_action=None, version="vectorised"):
        cdef Pool mem = Pool()
        cdef Py_ssize_t length, len_x, i, j, n, i_bcs, i_save
        cdef double[::1] getarr = np.zeros(2, np.float)
        # array to store the wave speeds
        cdef double[::1] c_i
        cdef CFL res_cfl
        cdef SOLVER res_solver
        length = len(self.mesh.vessels)

        cdef Py_ssize_t siz_i
        cdef vector[U] v_U
        cdef vector[vector[Solution]] vSol
        v_U.resize(length)
        self.initialise_solution_vector(v_U, length)

        import time
        t0 = time.perf_counter()  # CPU time measurement

        # --- Valid indices for space and time mesh ---
        self._It = np.arange(0, self._Nt + 1, dtype=np.int64)

        self._dt = self._t[1] - self._t[0]

        print("Solver set to dt=%0.9f" % self._dt)

        cdef vector[cParabolicSolObj] parU
        parU.resize(length)
        for i in range(length):
            N_n = self.mesh.vessels[i].r0.shape[0]
            parU[i]._size = N_n
            parU[i].x_hyper = <double**>self.mem.alloc(2, sizeof(double*))
            for j in range(2):
                parU[i].x_hyper[j] = <double*>self.mem.alloc(N_n, sizeof(double))
            parU[i].a = <double*>self.mem.alloc(N_n-3, sizeof(double))
            parU[i].b = <double*>self.mem.alloc(N_n-2, sizeof(double))
            parU[i].c = <double*>self.mem.alloc(N_n-3, sizeof(double))
            parU[i].d = <double*>self.mem.alloc(N_n, sizeof(double))
            parU[i].a_ph = <double*>self.mem.alloc(N_n-2, sizeof(double))
            parU[i].a_mh = <double*>self.mem.alloc(N_n-2, sizeof(double))

        if user_action is not None:
            i_save = 0
            vSol.push_back(vector[Solution]())
            vSol[i_save].resize(length)

        # check with nogil here
        for i in range(length):
            len_x = self.mesh.vessels[i].x.shape[0]
            if user_action is not None:
                vSol[i_save][i].n = 0
                vSol[i_save][i].size = 3*len_x
                vSol[i_save][i].u_store = <double*>mem.alloc(3*len_x, sizeof(double)) 
            for j in range(len_x):
                self._bcs.I(self.mesh.vessels[i].x[j], j, i, getarr)
                v_U[i].u_n[0][j] = getarr[0]
                v_U[i].u_n[1][j] = getarr[1]

                # check cfl condition
                res_cfl = self.cfl_i(getarr, v_U[i].dx, self._dt, j, i)
                if res_cfl == CFL.CFL_STATUS_ERROR:
                    return SOLVER.SOLVER_STATUS_ERROR

                if user_action is not None:
                    p = self._pdes.pressure_i(getarr[0], j, vessel_index=i)

                    vSol[i_save][i].u_store[j] = v_U[i].u_n[0][j]
                    vSol[i_save][i].u_store[j + len_x] = v_U[i].u_n[1][j]
                    vSol[i_save][i].u_store[j + 2*len_x] = p

#         # TIME LOOP
        cdef ITYPE_t index, siz_U0, k, siz_UL, siz_UBif, siz_UConj, par, d1, d2
        cdef double[::1] u_nmp1, x_n, u_nmp1_conj, x_n_conj
        cdef double* p_visco
        cdef Py_ssize_t ic
        u_nmp1 = np.zeros(6, np.float)
        x_n = np.zeros(6, np.float)
        u_nmp1_conj = np.zeros(4, np.float)
        x_n_conj = np.zeros(4, np.float)
        cdef ITYPE_t last_index, first_index, iter_count
        cdef ITYPE_t[::1] bif_list = np.zeros(3, dtype=np.int64)
        cdef ITYPE_t[::1] conj_list = np.zeros(2, dtype=np.int64)
        siz_U0 = self.U0_arr.shape[0]
        siz_UL = self.UL_arr.shape[0]
        siz_UBif = self.UBif_arr.shape[0]
        siz_UConj = self.UConj_arr.shape[0]

        #range [1 , self._Nt + 1]
        for n in range(1, self._Nt + 1):
            # this needs to be parallel
            # with nogil, parallel(num_threads=8):
            for i in prange(length, nogil=True):
            # for i in range(length):
                self.advance(v_U[i], self._dt, v_U[i].dx, i)

#             # Insert boundary conditions
            for index in range(siz_U0):
                k = self.U0_arr[index]
                self._bcs.U_0(v_U[k].u_n, v_U[k].x_size, self._t[n], v_U[k].dx, self._dt, k, v_U[k].u)

            # with nogil, parallel(num_threads=8):
            for index in prange(siz_UL, nogil=True):
                self._bcs.U_L(v_U[self.UL_arr[index]].u_n, v_U[self.UL_arr[index]].x_size, self._t[n], v_U[self.UL_arr[index]].dx,
                              self._dt, self.UL_arr[index], v_U[self.UL_arr[index]].u)

            #parallel siz_UBif
            for index in range(siz_UBif):
            # for index in prange(siz_UBif, nogil=True):
                par = self.UBif_arr[index, 0]
                d1 = self.UBif_arr[index, 1]
                d2 = self.UBif_arr[index, 2]
                bif_list[0] = par
                bif_list[1] = d1
                bif_list[2] = d2
                u_nmp1[0], u_nmp1[1], u_nmp1[2], u_nmp1[3], u_nmp1[4], u_nmp1[5] = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
                x_n[0], x_n[1], x_n[2], x_n[3], x_n[4], x_n[5] = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
                last_index = v_U[par].x_size
                first_index = 0
                u_nmp1[0] = v_U[par].u[0][last_index - 2]
                u_nmp1[1] = v_U[par].u[1][last_index - 2]
                # print v_U[par].u[1, last_index - 2]
                x_n[0] = v_U[par].u_n[0][last_index - 1]
                x_n[1] = v_U[par].u_n[1][last_index - 1]
                u_nmp1[2] = v_U[d1].u[0][first_index + 1]
                u_nmp1[3] = v_U[d1].u[1][first_index + 1]
                x_n[2] = v_U[d1].u_n[0][first_index]
                x_n[3] = v_U[d1].u_n[1][first_index]
                u_nmp1[4] = v_U[d2].u[0][first_index + 1]
                u_nmp1[5] = v_U[d2].u[1][first_index + 1]
                x_n[4] = v_U[d2].u_n[0][first_index]
                x_n[5] = v_U[d2].u_n[1][first_index]


                # return SOLVER.SOLVER_STATUS_ERROR
                iter_count = Newton_Raphson_Sys_bif(self._bcs, x_n, u_nmp1, self._dt,
                                                    bif_list, 1.0e-9, 100)

                # print iter_count
                v_U[par].u[0][last_index - 1] = x_n[0]
                v_U[par].u[1][last_index - 1] = x_n[1]
                v_U[d1].u[0][0] = x_n[2]
                v_U[d1].u[1][0] = x_n[3]
                v_U[d2].u[0][0] = x_n[4]
                v_U[d2].u[1][0] = x_n[5]

            for index in range(siz_UConj):
                d1 = self.UConj_arr[index, 0]
                d2 = self.UConj_arr[index, 1]
                conj_list[0] = d1
                conj_list[1] = d2
                u_nmp1_conj[0], u_nmp1_conj[1], u_nmp1_conj[2], u_nmp1_conj[3] = 0.0, 0.0, 0.0, 0.0
                x_n_conj[0], x_n_conj[1], x_n_conj[2], x_n_conj[3] = 0.0, 0.0, 0.0, 0.0
                last_index = v_U[d1].x_size
                first_index = 0
                u_nmp1_conj[0] = v_U[d1].u[0][last_index - 2]
                u_nmp1_conj[1] = v_U[d1].u[1][last_index - 2]
                x_n_conj[0] = v_U[d1].u_n[0][last_index - 1]
                x_n_conj[1] = v_U[d1].u_n[1][last_index - 1]
                u_nmp1_conj[2] = v_U[d2].u[0][first_index + 1]
                u_nmp1_conj[3] = v_U[d2].u[1][first_index + 1]
                x_n_conj[2] = v_U[d2].u_n[0][first_index]
                x_n_conj[3] = v_U[d2].u_n[1][first_index]

                iter_count = Newton_Raphson_Sys_conj(self._bcs, x_n_conj, u_nmp1_conj, self._dt,
                                                     conj_list, 1.0e-9, 100)

                v_U[d1].u[0][last_index - 1] = x_n_conj[0]
                v_U[d1].u[1][last_index - 1] = x_n_conj[1]
                v_U[d2].u[0][0] = x_n_conj[2]
                v_U[d2].u[1][0] = x_n_conj[3]


            if (user_action is not None) and (n % user_action.skip_frame == 0):
                i_save += 1
                vSol.push_back(vector[Solution]())
                vSol[i_save].resize(length)

            for i in range(length):
                copy2darrays(v_U[i].u, parU[i].x_hyper, 2, parU[i]._size)
                self._advance_Godunov_splitting(parU[i], v_U[i].u, v_U[i].x_size, n, self._dt, self.mesh.vessels[i].dx, i, 1.0)
                if (user_action is not None) and (n % user_action.skip_frame == 0):
                    # # check cfl condition
                    res_cfl = self.cfl_condition(parU[i].x_hyper, parU[i]._size, v_U[i].dx, self._dt, i)
                    if res_cfl == CFL.CFL_STATUS_ERROR:
                        print("Solver failed in vessel %d in time increment %d  spatial %d and lenght=%d" % (i, n, j, v_U[i].x_size))
                        return SOLVER.SOLVER_STATUS_ERROR
                    else:
                        vSol[i_save][i].n = n
                        vSol[i_save][i].size = 3*v_U[i].x_size
                        vSol[i_save][i].u_store = <double*>mem.alloc(vSol[i_save][i].size, sizeof(double)) 
                        p_visco = <double*>mem.alloc(v_U[i].x_size, sizeof(double))
                        self._pdes.pressure_visco(v_U[i].u, v_U[i].x_size, i, p_visco)

                        for ic in range(v_U[i].x_size):
                            vSol[i_save][i].u_store[ic] = v_U[i].u[0][ic]
                            vSol[i_save][i].u_store[ic + v_U[i].x_size] = v_U[i].u[1][ic]
                            vSol[i_save][i].u_store[ic + 2*v_U[i].x_size] = p_visco[ic]

                for j in range(v_U[i].x_size):
                    v_U[i].u_n[0][j] = v_U[i].u[0][j]
                    v_U[i].u_n[1][j] = v_U[i].u[1][j]


            printf("\r%6.2f%%", (100.0 * self._t[n] / self._T))
        
        cdef Py_ssize_t size_solution = vSol.size()
        cdef Py_ssize_t id, jd
        cdef double[:, ::1] u_store
        if user_action is not None:

            print("Size of solution is:", size_solution)
            print("Writing solution to file...")
            
            for id in range(size_solution):
                for jd in range(length):
                    u_store = <double[:3, :vSol[id][jd].size/3]>(vSol[id][jd].u_store)
                    user_action(u_store, self.mesh.vessels[jd].x, self._t, vSol[id][jd].n, PRINT_STATUS, WRITE_STATUS, jd)
            print("Writing solution to file completed!")

        return SOLVER.SOLVER_STATUS_OK


cdef void copy2darrays(double** x, double** out, Py_ssize_t rows, Py_ssize_t clms):
    cdef Py_ssize_t i, j
    for i in range(rows):
        for j in range(clms):
            out[i][j] = x[i][j]

@cython.boundscheck(False)
cpdef double norm2(double[::1] x):
    cdef Py_ssize_t siz, i
    cdef double sums
    sums = 0.
    siz = x.shape[0]
    for i in range(siz):
        sums = sums + x[i]*x[i]
    return sqrt(sums)

@cython.boundscheck(False)
cpdef double dotVec(double[::1] a, double[::1] b):
    cdef Py_ssize_t siz, i
    cdef double sums
    siz = a.shape[0]
    sums = 0
    for i in range(siz):
        sums = sums + a[i]*b[i]
    return sums

@cython.boundscheck(False)
cdef void swaprows(double[:, ::1] a, ITYPE_t i, ITYPE_t j):
    cdef ITYPE_t clms = a.shape[1]
    cdef ITYPE_t k
    cdef double temp

    for k in range(clms):
        temp = a[i, k]
        a[i, k] = a[j, k]
        a[j, k] = temp

@cython.boundscheck(False)
cdef void swaprowsVd(double[::1] a, ITYPE_t i, ITYPE_t j):
    cdef double temp

    temp = a[i]
    a[i] = a[j]
    a[j] = temp

@cython.boundscheck(False)
cdef void swaprowsVi(ITYPE_t[::1] a, ITYPE_t i, ITYPE_t j):
    cdef ITYPE_t temp

    temp = a[i]
    a[i] = a[j]
    a[j] = temp

@cython.boundscheck(False)
@cython.cdivision(True)
cdef int cLUdecomp(double[:,::1] a, ITYPE_t[::1] seq, double tol=1.0e-09):
    """LU Doolittle's decomposition"""
    cdef ITYPE_t n, i, k, j, max_i, iterat
    cdef double lam, max_value
    cdef double[::1] s

    n = a.shape[0]
    s = np.zeros(n, dtype=np.float)

    for i in range(n):
        max_value = 0.0
        seq[i] = i
        for j in range(a.shape[1]):
            if (fabs(a[i, j]) > max_value):
                max_value = fabs(a[i, j])
        s[i] = max_value

    for k in range(0, n-1):
        max_value = 0.0
        max_i = 0
        iterat = 0
        for j in range(k, n):
            if (abs(a[j,k])/s[j] > max_value):
                max_value = abs(a[j,k])/s[j]
                max_i = iterat
            iterat= iterat + 1
        p = max_i + k

        if abs(a[p, k]) < tol:
            return -1
        if p != k:
            swaprowsVd(s, k, p)
            swaprows(a, k, p)
            swaprowsVi(seq, k, p)

        for i in range(k+1, n):
            if a[i, k] != 0.0:
                lam = a[i, k] / a[k, k]
                for j in range(k+1, n):
                    a[i, j] = a[i, j] - lam*a[k, j]
                a[i, k] = lam
    return 0

@cython.cdivision(True)
@cython.boundscheck(False)
cdef ITYPE_t cLUsolve(double[:, ::1] a, double[::1] b, ITYPE_t[::1] seq, double[::1] out):
    cdef Py_ssize_t n, siz, i, j, k
    cdef double sums

    cLUdecomp(a, seq)

    n = a.shape[0]
    for i in range(out.shape[0]):
        out[i] = b[seq[i]]

    for k in range(1, n):
        sums = 0.0
        for i in range(0, k):
            sums = sums + a[k, i]*out[i]
        out[k] = out[k] - sums

    for j in range(n):
        k = n - 1 - j
        sums = 0.0
        for i in range(k+1, n):
            sums = sums + a[k, i]*out[i]
        out[k] = (out[k] - sums)/a[k, k]
    return 0

@cython.initializedcheck(False)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef ITYPE_t Newton_Raphson_Sys_bif(cBCs bcs, double[::1] x, double[::1] u, double dt, ITYPE_t[::1] vessel_indices,
                                double eps, int N=100):
    cdef Py_ssize_t n, siz, i, j
    cdef ITYPE_t iter_count
    cdef double norm

    cdef double[::1] R_value, delta
    cdef double[:, ::1] J_value
    cdef ITYPE_t[::1] seq

    n = x.shape[0]

    R_value = np.zeros(n, np.float)
    delta = np.zeros(n, np.float)
    J_value = np.zeros((n, n), dtype=np.float)
    seq = np.zeros(n, dtype=np.int64)

    bcs.bifurcation_R(x, u, dt, vessel_indices, R_value)
    norm = norm2(R_value)

    # we need to solve to negative R value
    for i in range(n):
        R_value[i] = -R_value[i]

    iter_count = 0
    while fabs(norm) > eps and iter_count < N:
        bcs.bifurcation_Jr(x, vessel_indices, J_value)
        cLUsolve(J_value, R_value, seq, delta)
        for i in range(n):
            x[i] = x[i] + delta[i]
        bcs.bifurcation_R(x, u, dt, vessel_indices, R_value)
        norm = norm2(R_value)
        for i in range(n):
            R_value[i] = -R_value[i]
        iter_count += 1
    if fabs(norm) > eps:
        iter_count = -1
    return iter_count

@cython.initializedcheck(False)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef int Newton_Raphson_Sys_conj(cBCs bcs, double[::1] x, double[::1] u, double dt, ITYPE_t[::1] vessel_indices,
                                 double eps, int N=100):
    cdef Py_ssize_t n, siz, i
    cdef ITYPE_t iter_count
    cdef double norm

    cdef double[::1] R_value, delta
    cdef double[:, ::1] J_value
    cdef ITYPE_t[::1] seq
    n = x.shape[0]

    R_value = np.zeros(n, np.float)
    delta = np.zeros(n, np.float)
    J_value = np.zeros((n, n), dtype=np.float)
    seq = np.zeros(n, dtype=np.int64)

    bcs.conjuction_R(x, u, dt, vessel_indices, R_value)

    norm = norm2(R_value)

    #we need to solve to negative R value
    for i in range(n):
        R_value[i] = -R_value[i]
    iter_count = 0
    while fabs(norm) > eps and iter_count < N:
        bcs.conjuction_Jr(x, vessel_indices, J_value)
        cLUsolve(J_value, R_value, seq, delta)
        for i in range(n):
            x[i] = x[i] + delta[i]
        bcs.conjuction_R(x, u, dt, vessel_indices, R_value)
        norm = norm2(R_value)
        for i in range(n):
            R_value[i] = -R_value[i]
        iter_count += 1
    if fabs(norm) > eps:
        iter_count = -1
    return iter_count

@cython.cdivision(True)
cpdef double linear_extrapolation(double x, double x1, double x2, double y1, double y2):
    return y1 + (x - x1)*((y2 - y1)/(x2-x1))