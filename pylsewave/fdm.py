"""
pylsewave module containing classes for finite difference numerical computation.
E.g. MacCormack, Lax-Wedroff, etc.

Parts of this module have been recycled from `Finite difference computing with
partial differential equations <https://github.com/hplgit/fdm-book>`_
"""
from __future__ import division

import numpy as np
import os
from timeit import default_timer as timer
# from progress import progress_timer
from pylsewave.mesh import VesselNetwork
from pylsewave.nonlinearsolvers import Newton_system_conj_points
import sys
from sys import stdout
from pylsewave.pwutils import CV_fun, h_walls
from pylsewave.cynum import pytdma
from pylsewave.pdes import PDEm
from pylsewave.bcs import BCs
# from scipy.optimize import fsolve


__author__ = "Georgios E. Ragkousis"
__acknowledgement__ = "Parts have been recycled from Hans Peter Langtangen"

WRITE_STATUS = True
PRINT_STATUS = False


class FDMSolver(object):
    """
    Base class for finite-difference computing schemes.
    """
    def __init__(self, inbcs):
        """
        :param pylsewave.bcs.BCs inbcs: class with the bcs definition
        :raises TypeError: if the in pdes is not type pylsewave.pdes.PDEm
        """
        if isinstance(inbcs, BCs) is not True:
            raise TypeError("The input PDEs should be of type pylsewave.bcs.BCs!")
        self._bcs = inbcs
        self._pdes = inbcs.pdes
        self.mesh = inbcs.mesh
        self._no_cycles = None
        self._T = None
        self._dt = None
        self._Nt = None
        self._t = None
        self.U0_arr = None
        self.UL_arr = None
        self.UBif_arr = None
        self.UConj_arr = None

    @property
    def dt(self):
        """

        :getter: returns the dt of the solver
        :type: float
        """
        return self._dt

    def set_T(self, dt, T, no_cycles):
        """
        Method to set the total period along with the time stepping of the solver.
        :param float dt: time step
        :param float T: total period
        :param no_cycles: number of the cycles (if periodic)
        :return: None
        """
        self._no_cycles = no_cycles
        self._T = T
        self._dt = dt
        self._Nt = int(round(self._T / self._dt))
        self._t = np.linspace(0.0, self._Nt * self._dt, self._Nt + 1)
        self._pdes.set_boundary_layer_th(T, no_cycles)

    def set_BC(self, U0_array,
               UL_array, UBif_array=None, UConj_array=None):
        """
        Method to set the BCs/connectivity of the network/mesh

        :param numpy.ndarray U0_array: array with vessel Ids that inlet bcs will be prescribed
        :param numpy.ndarray UL_array: array with vessel Ids that outlet bcs will be prescribed
        :param numpy.ndarray UBif_array: array with vessel Ids forming bifurcations
        :param numpy.ndarray UConj_array: array with vessel Ids forming conjunctions
        :return: None
        """
        self.U0_arr = U0_array
        self.UL_arr = UL_array
        self.UBif_arr = UBif_array
        self.UConj_arr = UConj_array

    def cfl_condition(self, u, dx, dt, vessel_index):
        """
        Method to check the pre-defined CFL condition.

        :param numpy.ndarray u: array with the solution u
        :param dx: spatial discretisation of the vessel
        :param dt: time step
        :param vessel_index: unique vessel Id
        :return: bool
        """
        _A, _q = u
#         print u
        _theta = dt / dx
        _c = self._pdes.compute_c(_A, None, vessel_index)
#         print _c
        _u = _q / _A
#         _W = [_u + _c, _u - _c]
        _W = _u + _c
        _W_plus = np.absolute(_W)
        _fraction = 1. / _W_plus
        comparison_elem = np.amin(_fraction)
        if (_theta < comparison_elem) == True:
            return True
        else:
            return False

    def solve(self, casename, user_action=None,
              version='vectorized'):
        """
        Method to solve the hyperbolic part of the PDE system:

        .. math::

            U_t + F(U)_x = S(U)

        :param str casename: the name of the simulation case
        :param user_action: a callback function/class to store the solution
        :param str version: scalar or vectorised computations
        :return: int
        """
        # append spatial mesh to list
        length = len(self.mesh.vessels)
        self._dt = self._t[1] - self._t[0]

        # --- Make hash of all input data ---
        # clever way of checking whether the analysis has run ag
        # ain
        import hashlib
        import inspect
        data = (str(self._dt) + '_' +
                str(self._T) + '_').encode('utf-8')
        hashed_input = casename + hashlib.sha1(data).hexdigest()
        if os.path.isfile('.' + hashed_input + '_archive.npz'):
            # Simulation is already run
            return -1, hashed_input

        # --- Allocate memomry for solutions ---
        u_store = [np.zeros((3, len(self.mesh.vessels[i].x))) for i in range(length)]
        u = [np.zeros((2, len(self.mesh.vessels[i].x))) for i in range(length)]
        u_n = [np.zeros((2, len(self.mesh.vessels[i].x))) for i in range(length)]

        # --- Allocate memomry for solutions ---
        u_store = [np.zeros((3, len(self.mesh.vessels[i].x))) for i in range(length)]
        u = [np.zeros((2, len(self.mesh.vessels[i].x))) for i in range(length)]
        u_n = [np.zeros((2, len(self.mesh.vessels[i].x))) for i in range(length)]

        for i in range(length):
            u_n[i] = self._bcs.I(self.mesh.vessels[i].x, i)

            if user_action is not None:
                p = self._pdes.pressure(u_n[i][0, :], vessel_index=i)
                u_store[i][0:-1, :] = u_n[i].copy()
                u_store[i][-1, :] = p.copy()
                user_action(u_store[i], self.mesh.vessels[i].x, self._t, 0, PRINT_STATUS, WRITE_STATUS, i)

        print("Solver set to dt=%0.9f" % self._dt)

        self._It = range(0, self._Nt + 1)
        n = 1
        while n < self._Nt + 1:
#         for n in self._It[1:]:
            for i in range(length):
                self._advance(u[i], u_n[i], n, self._dt, self.mesh.vessels[i].dx, i)

            #BCS
            for k in self.U0_arr:
                self._bcs.U_0(u_n[k], self._t[n], self.mesh.vessels[k].dx, self._dt, k, u[k])

            for k in self.UL_arr:
                self._bcs.U_L(u_n[k], self._t[n], self.mesh.vessels[k].dx, self._dt, k, u[k])

            if self.UBif_arr is not None:
                for i in range(len(self.UBif_arr)):
                    p, d1, d2 = self.UBif_arr[i]
                    u_nmp1 = np.zeros(6)
                    x_n = np.zeros(6)
                    u_nmp1[0:2] = u[p][:, -2]
                    x_n[0:2] = u_n[p][:, -1]
                    u_nmp1[2:4] = u[d1][:, 1]
                    x_n[2:4] = u_n[d1][:, 0]
                    u_nmp1[4:6] = u[d2][:, 1]
                    x_n[4:6] = u_n[d2][:, 0]
                    x, iters = Newton_system_conj_points(self._bcs.bifurcation_R, self._bcs.bifurcation_Jr,
                                                         x_n, u_nmp1, self._dt,
                                                         [p, d1, d2], eps=1.0e-9,
                                                         N=100)
                    # or you can use fsolve from scipy
                    # x = fsolve(self._bcs.bifurcation_R, x_n, args=(u_nmp1, self._dt, [p, d1, d2]),
                    #            maxfev=100, xtol=1.0e-9)
                    u[p][:, -1] = x[0:2]
                    u[d1][:, 0] = x[2:4]
                    u[d2][:, 0] = x[4:6]

            if self.UConj_arr is not None:
                for i in range(len(self.UConj_arr)):
                    d1, d2 = self.UConj_arr[i]
                    u_nmp1 = np.zeros(4)
                    x_n = np.zeros(4)
                    u_nmp1[0:2] = u[d1][:, -2]
                    x_n[0:2] = u_n[d1][:, -1]
                    u_nmp1[2:4] = u[d2][:, 1]
                    x_n[2:4] = u_n[d2][:, 0]
                    x, iters = Newton_system_conj_points(self._bcs.conjuction_R, self._bcs.conjuction_Jr,
                                                         x_n, u_nmp1, self._dt,
                                                         [d1, d2], eps=1.0e-9,
                                                         N=100)
                    # x = fsolve(self._bcs.conjuction_R, x_n, args=(u_nmp1, self._dt, [d1, d2]),
                    #            maxfev=100, xtol=1.0e-9)
                    u[d1][:, -1] = x[0:2]
                    u[d2][:, 0] = x[2:4]

            for i in range(length):

                if (user_action is not None) and (n % user_action.skip_frame == 0):
                    # check CFL condition
                    res = self.cfl_condition(u[i], self.mesh.vessels[i].dx, self._dt, i)
                    if res is True:
                        # store solution
                        p = self._pdes.pressure(u[i][0, :], None, i)
                        u_store[i][0:-1, :] = u[i].copy()
                        u_store[i][-1, :] = p.copy()
                        user_action(u_store[i], self.mesh.vessels[i].x, self._t, n, PRINT_STATUS, WRITE_STATUS, i)
                    else:
                        print("Solver failed in vessel %d in time increment %d " % (i, n))
                        raise ValueError("The CFL condition hasn't been satisfied!\n"
                                         "Reduce dt")
                # if splitting_scheme == 'Godunov':
                #     self._Godunov_splitting(u[i], n, self._dt, self._dt, self.mesh.vessels[i].dx, i)
                # here we transfer the solution back to u_n arrays
                u_n[i] = u[i].copy()


            # progress = (100.0 * self._t[n] / self._T)
            sys.stdout.write('\r')
            sys.stdout.write("[%-50s] %6.2f%%" % ('='*int(round((100.0 * self._t[n] / self._T)/2.0)),
                                                  (100.0 * self._t[n] / self._T)))
            # stdout.write("\r%6.2f %%" % (100.0 * self._t[n] / self._T))
            stdout.flush()

            n += 1


        print("The time increment dt was: %0.9f" % self._dt)

        return STATUS_OK

    def _advance(self, u, u_n, n, dt, dx, vessel_index):
        """
        Method to advance the solution at every time step according to the FDM scheme.

        :param numpy.ndarray u: the solution vector at time n+1
        :param numpy.ndarray u_n: the solution vector at time n
        :param int n: the number of time step
        :param float dt: the time step
        :param float dx: the spatial discretisation
        :param vessel_index:  the unique vessel index
        :return: None
        :raises NotImplementedError: use a sub-class and override this method
        """
        raise NotImplementedError


class BloodWaveLaxWendroff(FDMSolver):
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
    def _advance(self, u, u_n, n, dt, dx, vessel_index):

        theta = dt/dx
        dt2 = dt*0.5

        u_nh_mh = np.zeros((2, u.shape[1]))  # predictor
        u_nh_ph = np.zeros((2, u.shape[1]))  # predictor

        F_ = self._pdes.flux(u_n, 0, None, vessel_index)
        S_ = self._pdes.source(u_n, 0, None, vessel_index)

        # U[n+1/2, i-1/2]
        u_nh_mh[:, 1:-1] = ((u_n[:, 1:-1] + u_n[:, 0:-2]) / 2. -
                               0.5 * (theta) * (F_[:, 1:-1] - F_[:, 0:-2]) +
                               0.5 * dt2 * (S_[:, 1:-1] + S_[:, 0:-2]))
        # U[n+1/2, i+1/2]
        u_nh_ph[:, 1:-1] = ((u_n[:, 2:] + u_n[:, 1:-1]) / 2. -
                               0.5 * (theta) * (F_[:, 2:] - F_[:, 1:-1]) +
                               0.5 * dt2 * (S_[:, 2:] + S_[:, 1:-1]))
        # this is a test, because we need to fill the edges with ghost
        u_nh_mh[:, 0] = u_n[:, 0]
        u_nh_mh[:, -1] = u_n[:, -1]
        u_nh_ph[:, 0] = u_n[:, 0]
        u_nh_ph[:, -1] = u_n[:, -1]

        F_nh_mh = self._pdes.flux(u_nh_mh, -1, None, vessel_index)
        S_nh_mh = self._pdes.source(u_nh_mh, -1, None, vessel_index)
        F_nh_ph = self._pdes.flux(u_nh_ph, 1, None, vessel_index)
        S_nh_ph = self._pdes.source(u_nh_ph, 1, None, vessel_index)

        u[:, 1:-1] = (u_n[:, 1:-1] - theta * (F_nh_ph[:, 1:-1] - F_nh_mh[:, 1:-1]) +
                      dt2*(S_nh_ph[:, 1:-1] + S_nh_mh[:, 1:-1]))


class BloodWaveMacCormack(FDMSolver):
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
    def _advance(self, u, u_n, n, dt, dx, vessel_index):

        theta = (dt / dx)
        dt2 = dt*0.5

        u_star = np.zeros((2, u.shape[1]))  # predictor

        F_ = self._pdes.flux(u_n, 0, None, vessel_index)
        S_ = self._pdes.source(u_n, 0, None, vessel_index)
        # U*
        u_star[:, :-1] = (u_n[:, :-1] - theta*(F_[:, 1:] - F_[:, :-1]) +
                          dt*S_[:, :-1])
        # fill the arrays
        # this will fill the last element in array
        u_star[:, -1] = (u_n[:, -1] - theta*(F_[:, -2] - F_[:, -1]) + dt*S_[:, -1])

        F_star = self._pdes.flux(u_star, 0, None, vessel_index)
        S_star = self._pdes.source(u_star, 0, None, vessel_index)

        u[:, 1:] = (0.5 *(u_n[:, 1:] + u_star[:, 1:] - theta*(F_star[:, 1:] -
                    F_star[:, :-1]) + dt*S_star[:, 1:]))


class BloodWaveMacCormackGodunov(BloodWaveMacCormack):

    def solve(self, casename, user_action=None,
              version='vectorized'):
        """
        Solve U_t + F_x = S
        """
        # append spatial mesh to list

        length = len(self.mesh.vessels)
        self._dt = self._t[1] - self._t[0]
        # a small Qt app with progressbar
        #         self._pdes.dt = dt
        # if progress_bar is not False:
        #     from pwqtutils import pylseWaveProgressBar
        #     from PyQt5.QtWidgets import (QApplication)
        #
        #     app = QApplication(sys.argv)
        #     pbar = pylseWaveProgressBar(dt=self._dt, T=self._T)
        #     pbar.show()
        #     pbar.progressLabel.setText("Calculating the Godunov splitting matrices...")
        #     app.processEvents()
        # --- Make hash of all input data ---
        # clever way of checking whether the analysis has run ag
        # ain
        import hashlib
        import inspect
        data = (str(self._dt) + '_' +
                str(self._T) + '_').encode('utf-8')
        hashed_input = casename + hashlib.sha1(data).hexdigest()
        if os.path.isfile('.' + hashed_input + '_archive.npz'):
            # Simulation is already run
            return -1, hashed_input

        # --- Allocate memomry for solutions ---
        u_store = [np.zeros((3, len(self.mesh.vessels[i].x))) for i in range(length)]
        u = [np.zeros((2, len(self.mesh.vessels[i].x))) for i in range(length)]
        u_n = [np.zeros((2, len(self.mesh.vessels[i].x))) for i in range(length)]

        # --- Allocate memomry for solutions ---
        u_store = [np.zeros((3, len(self.mesh.vessels[i].x))) for i in range(length)]
        u = [np.zeros((2, len(self.mesh.vessels[i].x))) for i in range(length)]
        u_n = [np.zeros((2, len(self.mesh.vessels[i].x))) for i in range(length)]

        # calculate all the vectors (lower, diagonal, upper)
        # self._set_vectors_Godunov_splitting(theta=1.0)

        for i in range(length):
            u_n[i] = self._bcs.I(self.mesh.vessels[i].x, i)

            if user_action is not None:
                p = self._pdes.pressure(u_n[i][0, :], u_n[i][1, :], vessel_index=i)
                u_store[i][0:-1, :] = u_n[i].copy()
                u_store[i][-1, :] = p.copy()
                user_action(u_store[i], self.mesh.vessels[i].x, self._t, 0, PRINT_STATUS, WRITE_STATUS, i)

        print("Solver set to dt=%0.9f" % self._dt)

        self._It = range(0, self._Nt + 1)
        # if progress_bar is not False:
        #     pbar.progressLabel.setText("Solving pylse wave problem...")
        #     app.processEvents()

        n = 1
        while n < self._Nt + 1:

            #         for n in self._It[1:]:
            for i in range(length):
                self._advance(u[i], u_n[i], n, self._dt, self.mesh.vessels[i].dx, i)

            # BCS
            for k in self.U0_arr:
                self._bcs.U_0(u_n[k], self._t[n], self.mesh.vessels[k].dx, self._dt, k, u[k])

            for k in self.UL_arr:
                self._bcs.U_L(u_n[k], self._t[n], self.mesh.vessels[k].dx, self._dt, k, u[k])

            if self.UBif_arr is not None:
                for i in range(len(self.UBif_arr)):
                    p, d1, d2 = self.UBif_arr[i]
                    u_nmp1 = np.zeros(6)
                    x_n = np.zeros(6)
                    u_nmp1[0:2] = u[p][:, -2]
                    x_n[0:2] = u_n[p][:, -1]
                    u_nmp1[2:4] = u[d1][:, 1]
                    x_n[2:4] = u_n[d1][:, 0]
                    u_nmp1[4:6] = u[d2][:, 1]
                    x_n[4:6] = u_n[d2][:, 0]
                    x, iters = Newton_system_conj_points(self._bcs.bifurcation_R, self._bcs.bifurcation_Jr,
                                                         x_n, u_nmp1, self._dt,
                                                         [p, d1, d2], eps=1.0e-9,
                                                         N=100)
                    u[p][:, -1] = x[0:2]
                    u[d1][:, 0] = x[2:4]
                    u[d2][:, 0] = x[4:6]

            if self.UConj_arr is not None:
                for i in range(len(self.UConj_arr)):
                    d1, d2 = self.UConj_arr[i]
                    u_nmp1 = np.zeros(4)
                    x_n = np.zeros(4)
                    u_nmp1[0:2] = u[d1][:, -2]
                    x_n[0:2] = u_n[d1][:, -1]
                    u_nmp1[2:4] = u[d2][:, 1]
                    x_n[2:4] = u_n[d2][:, 0]
                    x, iters = Newton_system_conj_points(self._bcs.conjuction_R, self._bcs.conjuction_Jr,
                                                         x_n, u_nmp1, self._dt,
                                                         [d1, d2], eps=1.0e-9,
                                                         N=100)
                    u[d1][:, -1] = x[0:2]
                    u[d2][:, 0] = x[2:4]

            for i in range(length):
                x = u[i].copy()
                self._advance_Godunov_splitting(u[i], n, self._dt, self.mesh.vessels[i].dx, i)
                if (user_action is not None) and (n % user_action.skip_frame == 0):
                    # check CFL condition
                    res = self.cfl_condition(x, self.mesh.vessels[i].dx, self._dt, i)
                    if res is True:
                        # store solution
                        p = self._pdes.pressure(u[i][0, :], u[i][1, :], None, i)
                        u_store[i][0:-1, :] = u[i].copy()
                        u_store[i][-1, :] = p.copy()
                        user_action(u_store[i], self.mesh.vessels[i].x, self._t, n, PRINT_STATUS, WRITE_STATUS, i)
                    else:
                        print("Solver failed in vessel %d in time increment %d " % (i, n))
                        raise ValueError("The CFL condition hasn't been satisfied!\n"
                                         "Reduce dt")

                u_n[i] = u[i].copy()


            progress = (100.0 * self._t[n] / self._T)
            stdout.write('\r')
            stdout.write("[%-50s] %6.2f%%" % ('='*int(round(progress/2.0)), progress))

            # stdout.write("\r%6.2f %%" % progress)
            stdout.flush()
            # pbar.pbar.setFormat("%6.2f%%" % progress)
            # # print(ex.step)
            # pbar.pbar.setValue(progress)
            # app.processEvents()
            n += 1

        print("The time increment dt was: %0.9f" % self._dt)

        return STATUS_OK

    def _advance_Godunov_splitting(self, u, n, dt, dx, vessel_index, theta=1.0):
        w_th = h_walls(self.mesh.vessels[vessel_index].r0)
        C_v = self._pdes.cv_f_qua(u[0, :], w_th, vessel_index)
        N_n = u.shape[1]
        dx = self.mesh.vessels[vessel_index].dx
        F = self._dt / (dx * dx)
        N_n = self.mesh.vessels[vessel_index].r0.shape[0]
        dia = np.zeros(N_n - 2)
        lower = np.zeros(N_n)
        lw = np.zeros(N_n - 3)
        upper = np.zeros(N_n)
        upr = np.zeros(N_n - 3)

        a_ph = 0.5 * (C_v[1:-1] + C_v[2:])
        a_mh = 0.5 * (C_v[1:-1] + C_v[:-2])

        dia = (1 + (u[0,1:-1]/self._pdes._rho)*theta * F * (a_ph + a_mh))
        lw[:] = -(u[0,2:-1]/self._pdes._rho)*theta * F * a_mh[1:]
        upr[:] = -(u[0,1:-2]/self._pdes._rho)*theta * F * a_ph[:-1]

        # # ---- BCd Neuman ----- ##
        dia[0] = (1 + (u[0, 1]/self._pdes._rho)*theta * F * (a_ph[0] + a_mh[0])) - (u[0, 1]/self._pdes._rho)*theta * F * a_mh[0]
        # here the value for m, m-1
        dia[-1] = (1 + (u[0, -2]/self._pdes._rho)*theta * F * (a_ph[-1] + a_mh[-1])) - (u[0, -2]/self._pdes._rho)*theta * F * a_ph[-1]

        b = np.zeros(N_n)
        # C_v = self.C_v[vessel_index]
        # a_ph = 0.5 * (C_v[1:-1] + C_v[2:])
        # a_mh = 0.5 * (C_v[1:-1] + C_v[:-2])
        # F = dt / (dx * dx)
        gamma = (u[0,1:-1]/self._pdes._rho)*(1 - theta)

        b[1:-1] = (gamma * F * a_mh * u[1, :-2] + (1 - gamma * F * (a_mh + a_ph)) * u[1, 1:-1] +
                   gamma * F * a_ph * u[1, 2:])

        res = pytdma(lw, dia, upr, b[1:-1], u[1, 1:-1])

if __name__ == '__main__':
    raise NotImplementedError('Module is not idented for direct execution')