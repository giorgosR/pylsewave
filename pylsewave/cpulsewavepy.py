""" Module with python function wrappers for c++ dll (ctypes) """
__author__ = "Georgios E. Ragkousis"
from ctypes import (cdll, CFUNCTYPE, c_double,
                    c_int, POINTER, c_float, pointer,
                    c_size_t)
from numpy.ctypeslib import ndpointer
import numpy as np

pwpydll = cdll.LoadLibrary(r'./CLibs/build/Debug/cpulsewavepy.dll')

doube_p1d = ndpointer(dtype=np.double, ndim=1, flags="CONTIGUOUS")
doube_p2d = ndpointer(dtype=np.uintp, ndim=1, flags="CONTIGUOUS")
doube_p3d = ndpointer(dtype=np.uintp, ndim=3, flags="CONTIGUOUS")

test_fun = pwpydll.sum
test_fun.restype = c_int
test_fun.argtypes = [c_int, c_int]

compute_radius_f = pwpydll.compute_radius
compute_radius_f.restype = None
compute_radius_f.argtypes = [c_double, c_double, c_double,
                             ndpointer(dtype=np.double, ndim=1, flags="CONTIGUOUS"),
                             ndpointer(dtype=np.double, ndim=1, flags="CONTIGUOUS"),
                             c_int]

compute_radius_s_f = pwpydll.compute_radius_s
compute_radius_s_f.restype = None
compute_radius_s_f.argtypes = [c_double, c_double, c_double, c_double, POINTER(c_double)]

compute_area_f = pwpydll.compute_area
compute_area_f.restype = None
compute_area_f.argtypes = [c_double, c_double, c_double,
                             ndpointer(dtype=np.double, ndim=1, flags="CONTIGUOUS"),
                             ndpointer(dtype=np.double, ndim=1, flags="CONTIGUOUS"),
                             c_int]

solver_Lax_Wendroff = pwpydll.solve_Lax_Wendroff
solver_Lax_Wendroff.restype = None
solver_Lax_Wendroff.argtypes = [doube_p2d, doube_p2d, doube_p1d, doube_p1d,
                                c_int, c_int, c_int, doube_p2d]


def addpy(a, b):
    return test_fun(a, b)


def compute_radius(r_p, r_d, l, x_in, x_out):
    if hasattr(x_in, "__len__"):
        return compute_radius_f(r_p, r_d, l, x_in, x_out, len(x_in))
    else:
        return compute_radius_s_f(r_p, r_d, l, x_in, x_out)


def compute_area(r_p, r_d, l, x_in, x_out):
    return compute_area_f(r_p, r_d, l, x_in, x_out, len(x_in))


def lax_wendroff_solve(u, u_n, x, t):
    m1, m2 = u.shape
    n = len(t)
    u_out = np.zeros((n, m1, m2))
    upp = (u.__array_interface__['data'][0]
           + np.arange(u.shape[0]) * u.strides[0]).astype(np.uintp)
    unpp = (u_n.__array_interface__['data'][0]
            + np.arange(u_n.shape[0]) * u_n.strides[0]).astype(np.uintp)
    u_outpp = (u_out.__array_interface__['data'][0]
               + np.arange(u_out.shape[0])*u_out.strides[0]
               + u_out.strides[1]*np.arange(u_out.shape[1])).astype(np.uintp)


    return solver_Lax_Wendroff(upp, unpp, x, t, m1, m2, n, u_outpp)
