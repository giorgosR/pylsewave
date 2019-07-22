__author__ = "Georgios E. Ragkousis"
import numpy as np
from pwconsts import *
from mesh import VesselNetwork


def convert_data_periodic(xdata, ydata, cycles=2, plot=False):
    """
    Function to convert data to periodic.
    :param xdata:
    :param ydata:
    :param cycles:
    :param plot:
    :return:
    """
    if np.ndim(xdata) > 1 and np.ndim(ydata) > 1:
        raise ValueError("Input arrays should be one dimensional!")

    y_stack = np.vstack([ydata[0:-1] for i in range(cycles)])
    x_stack = np.vstack([xdata[0:-1] + xdata[-1]*i for i in range(cycles)])
    y_periodic = np.zeros(len(y_stack.flatten()) + 1)
    x_periodic = np.zeros(len(x_stack.flatten()) + 1)

    cycle_time = xdata[-1]

    y_periodic[0:-1] = y_stack.flatten()
    x_periodic[0:-1] = x_stack.flatten()
    y_periodic[-1] = y_periodic[0]
    x_periodic[-1] = cycles * cycle_time

    if plot is True:
        import matplotlib.pyplot as plt
        plt.plot(x_periodic, y_periodic)
        plt.show()

    return x_periodic, y_periodic

def compute_c(R0, k, rho):
    k1, k2, k3 = k
    return np.sqrt((2/(3.*rho))*(k2*np.exp(k3*R0) + k1))

def CV_fun( r0, w_th, rho, phi):
    _phi = phi
    w_th = h_walls(r0)
    C_v = (2./3)*((_phi * w_th)/(rho*r0))
    return C_v

def cv_fun( r0, w_th, a, phi):
    _phi = phi
    w_th = h_walls(r0)
    C_v = (2./3)*((_phi * w_th)/(a*r0))
    return C_v

def h_walls(r0):
    alpha = 0.2802
    beta = -5.053*0.1
    gamma = 0.1324
    d = -0.1114*0.1
    return r0*(alpha*np.exp(beta*r0) + gamma*np.exp(d*r0))

def calculate_min_tCFL(mesh, rho, min_no_nodes=10):
    if isinstance(mesh, VesselNetwork) is not True:
        raise TypeError("Mesh parameter should be of type pylsewave.mesh.VesselNetwork!")
    # check CFL and set dx accordingly
    siz_ves = len(mesh.vessels)
    compare_l_c0 = []
    for i in range(siz_ves):
        c_max = np.max(compute_c(mesh.vessels[i].r0, mesh.vessels[i].k, rho))
        A = np.pi * (mesh.vessels[i].r_prox * mesh.vessels[i].r_prox)
        compare_l_c0.append(mesh.vessels[i].length / c_max)

    min_value = min(compare_l_c0)
    index_min_value = np.argmin(compare_l_c0)
    print("The min length to wave speed radio has been computed to Vessel: '%s' " %
          mesh.vessels[index_min_value].name)

    # Nx_i = 1
    min_time = []
    for i in range(siz_ves):
        Nx_i = min_no_nodes * np.floor((mesh.vessels[i].length /
                                        compute_c(mesh.vessels[i].r_prox, mesh.vessels[i].k, rho)) /
                                       (min_value))
        dx_i = mesh.vessels[i].length / Nx_i
        mesh.vessels[i].dx = dx_i
        min_time.append(dx_i / np.max(compute_c(mesh.vessels[i].r0, mesh.vessels[i].k, rho)))

    return min(min_time)

def linear_extrapolation(x, x1, x2, y1, y2):
    return y1 + (x - x1)*((y2 - y1)/(x2-x1))

if __name__ == '__main__':
    raise NotImplementedError('Module is not idented for direct execution')