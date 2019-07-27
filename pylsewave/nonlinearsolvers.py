"""
pylsewave module with non-linear solvers (e.g. Newton-Raphson)
"""
from __future__ import division
import numpy as np
import sys

__author__ = "Georgios E. Ragkousis"
__notes__ = ("Some code has been recycled from the book: 'Programming for computations' "
             "by Svein Linge & Hans Peter Langtangen ")


def Newton(f, dfdx, x, eps, return_x_list=False):
    f_value = f(x)
    iteration_counter = 0
    if return_x_list:
        x_list = []

    while abs(f_value) > eps and iteration_counter < 1000:
        try:
            x = x - float(f_value)/dfdx(x)
        except ZeroDivisionError:
            print("Error! - derivative zero for x = ", x)
            sys.exit(1)     # Abort with error

        f_value = f(x)
        iteration_counter += 1
        if return_x_list:
            x_list.append(x)

    # Here, either a solution is found, or too many iterations
    if abs(f_value) > eps:
        iteration_counter = -1  # i.e., lack of convergence

    if return_x_list:
        return x_list, iteration_counter
    else:
        return x, iteration_counter


def Newton_system(F, J, x, eps, N=100):
    """
    Function for solving non linear system with Newton Raphson

    :param F: Vector F (residuals)
    :param J: Jacobian
    :param x: solution vector
    :param eps: tolerance of convergence
    :param N: iteration number
    :return: solution vector and number of iterations
    """
    F_value = F(x)
    F_norm = np.linalg.norm(F_value, ord=2) # l2 norm of vector
    iteration_counter = 0
    while abs(F_norm) > eps and iteration_counter < N:
        delta = np.linalg.solve(J(x), -F_value)
        x = x + delta
        F_value = F(x)
        F_norm = np.linalg.norm(F_value, ord=2)
        iteration_counter += 1
    # Here, either a solution is found, or too many iterations
    if abs(F_norm) > eps:
        iteration_counter = -1
    return x, iteration_counter


def Newton_system_conj_points(R, J, x, u, dt, vessel_indices, eps, N=100):
    """
    Function to compute the x vector with fluxes and areas at bifurcations
    :param F: vector residual function R_i
    :param J: vector Jacobian function dR_i / dx_j
    :param x: x_i vector of the initial guess
    :param u: u_{L-1/0+1} values to calculate the characteristics
    :param vessel_indices: indices for bifurcations
    :param eps: tolerance for Newton Raphson convergence
    :param N: number of iterations
    :return: x predicted values and number of iterations
    """

    # print("py x: ", x)
    # print("py u: ", u)
    R_value = R(x, u, dt, vessel_indices)
    # print("pyRvalue: ", R_value)
    R_norm = np.linalg.norm(R_value, ord=2) # l2 norm of vector
    # print("pyR_norm", R_norm)
    iteration_counter = 0
    while abs(R_norm) > eps and iteration_counter < N:
        # print("J(x) : ", J(x, vessel_indices))
        delta = np.linalg.solve(J(x, u, dt, vessel_indices), -R_value)
        # print("pyDelta1:", delta )
        x = x + delta
        R_value = R(x, u, dt, vessel_indices)
        # print("py New R_value", R_value)
        R_norm = np.linalg.norm(R_value, ord=2)
        # print("R_norm1: ", R_norm)
        iteration_counter += 1
    # Here, either a solution is found, or too many iterations
    if abs(R_norm) > eps:
        iteration_counter = -1
    return x, iteration_counter

if __name__ == '__main__':
    raise NotImplementedError('Module is not idented for direct execution')
