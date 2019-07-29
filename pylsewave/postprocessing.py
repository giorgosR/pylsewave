"""
Module with post-processing classes.
"""
from __future__ import division
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib import cm

__author__ = "Georgios E. Ragkousis"
__acknowledgement__ = "surface plotting implemented/recycled from vamPy (https://github.com/akdiem/vampy)"

class ExtractUVectors(object):
    """
    Class containing the results of a pylsewave run.
    """
    def __init__(self, filename):
        """
        Constructor

        :param filename: the filename with the compressed results.
        """
        try:
            self._fls = np.load(filename)
        except:
            raise ValueError("Filename cannot be unzipped!")
        self._t = self._fls["t"]
        self._x = None
        self._meshgrid_X = None
        self._meshgrid_T = None

    @property
    def t(self):
        return self._t

    @property
    def odbf(self):
        return self._fls

    @property
    def x(self):
        return self._x

    @x.setter
    def x(self, x):
        self._x = x

    @property
    def meshgrid_X(self):
        return self._meshgrid_X

    @meshgrid_X.setter
    def meshgrid_X(self, X):
        self._meshgrid_X = X

    @property
    def meshgrid_T(self):
        return self._meshgrid_T

    @meshgrid_T.setter
    def meshgrid_T(self, T):
        self._meshgrid_T = T

    def getUVector(self, vessel_no = None, cycle=None, no_points=None):
        """
        method to get the result fields from specific vessel at certain cycle.

        :param vessel_no: the vessel number
        :param cycle: cycle for extracting results
        :param no_points: no points for interpolation
        :return: A (area), Q (flow), P (pressure), U (velocity) fields
        """
        A, q, p = [], [], []

        if (vessel_no is None) and (cycle is not None):
            self.x = self.odbf["x"]
            for key in self.odbf.keys():
                if (key is not "t") and (key is not "x"):
                    A.append(self.odbf[key][0, :])
                    q.append(self.odbf[key][1, :])
                    p.append(self.odbf[key][2, :])

        elif (vessel_no is not None) and (cycle is not None):
            self.x = self.odbf["x_vessel_%03d" % vessel_no]
            for key in self.odbf.keys():
                if (key is not "t") and (key is not "x") and\
                     ("_vessel_%03d" % vessel_no in key) and\
                     ("x_vessel_%03d" % vessel_no not in key):
                    A.append(self.odbf[key][0, :])
                    q.append(self.odbf[key][1, :])
                    p.append(self.odbf[key][2, :])

        A, q, p = np.asarray(A), np.asarray(q), np.asarray(p)
        A_cycle = A[len(self.t) - (len(self.t) // cycle)::].T
        q_cycle = q[len(self.t) - (len(self.t) // cycle)::].T
        p_cycle = p[len(self.t) - (len(self.t) // cycle)::].T
        u_cycle = q_cycle / A_cycle

        f = interpolate.interp2d(self.t[len(self.t) - (len(self.t) // cycle)::],
                                 self.x, q_cycle, kind="cubic")
        g = interpolate.interp2d(self.t[len(self.t) - (len(self.t) // cycle)::],
                                 self.x, p_cycle, kind="cubic")
        k = interpolate.interp2d(self.t[len(self.t) - (len(self.t) // cycle)::],
                                 self.x, A_cycle, kind="cubic")
        l = interpolate.interp2d(self.t[len(self.t) - (len(self.t) // cycle)::],
                                 self.x, u_cycle, kind="cubic")

        if no_points:
            selected_t = np.linspace(self.t[len(self.t) - (len(self.t) // cycle)], self.t[-1], no_points)
            selected_x = np.linspace(0., self.x[-1], no_points)
            T, X = np.meshgrid(selected_t, selected_x)
            Q = f(selected_t, selected_x)
            P = g(selected_t, selected_x)
            A = k(selected_t, selected_x)
            U = l(selected_t, selected_x)

        else:
            x = np.linspace(0., self.x[-1], len(self.t[len(self.t) - (len(self.t) // cycle)::]))
            T, X = np.meshgrid(self.t[len(self.t) - (len(self.t) // cycle)::], x)
            Q = f(self.t[len(self.t) - (len(self.t) // cycle)::], x)
            P = g(self.t[len(self.t) - (len(self.t) // cycle)::], x)
            A = k(self.t[len(self.t) - (len(self.t) // cycle)::], x)
            U = l(self.t[len(self.t) - (len(self.t) // cycle)::], x)

        self.meshgrid_T = T
        self.meshgrid_X = X


        return A, Q, P, U


class VizVessel(object):
    """
    Class for visualising a certain vessel
    """
    def __init__(self, filename, id):
        """
        Constructor

        :param filename: filename with the results
        :param id: id of the vessel
        """
        self._id = id
        self._odbf = ExtractUVectors(filename)
        self._A, self._q, self._p, self._u = None, None, None, None
        self._T = None
        self._X = None

    def set_sampling_params(self, sampling_points, cycle):
        """
        Method to set the visualisation parameters.
        :param sampling_points: no of sampling points
        :param cycle: cycle number
        :return: None
        """
        self._A, self._q, self._p, self._u = self._odbf.getUVector(vessel_no=self._id,
                                                                   cycle=cycle,
                                                                   no_points=sampling_points)
        self._T = self._odbf.meshgrid_T
        self._X = self._odbf.meshgrid_X

    def set_scaling(self, scaling_array):
        """
        Method to scale fields
        :param scaling_array: scaling array
        :return: None
        """
        self._A = self._A*scaling_array[0]
        self._q = self._q*scaling_array[1]
        self._p = self._p*scaling_array[2]
        self._u = self._u*scaling_array[3]
        self._T = self._T*scaling_array[4]
        self._X = self._X*scaling_array[5]

    def plot_3Dfield(self, names_of_fields):
        """
        Method to plot a specific field
        :param names_of_fields: name of fields (Area, Flow, Pressure, Velocity, PU-LOOP)
        :return:
        """
        if "Area" in names_of_fields:
            # AREA
            fig3 = plt.figure(figsize=(8, 6), facecolor='w', edgecolor='k',
                              linewidth=2.0, frameon=True)
            ax = fig3.gca(projection='3d')
            surf = ax.plot_surface(self._T, self._X, self._A, rstride=4, cstride=4, alpha=0.2,
                                   cmap=cm.coolwarm, linewidth=0.2, antialiased=True)
            # ax.view_init(30, 45)
            # cset = ax.contour(T, X, A * rc ** 2, zdir='z', offset=0.0, cmap=cm.coolwarm, alpha=0.2)
            ax.set_xlabel('t')
            ax.set_ylabel('z')
            ax.set_zlabel('Area')
            ax.set_title('Area')
            # ax.set_xlim([min(t_scaled[len(t_scaled)-(len(t_scaled)/4)::])-1, max(t_scaled[len(t_scaled)-(len(t_scaled)/4)::])])
            # ax.set_ylim([min(x_scaled), 60])
            # fig3.colorbar(surf, shrink=0.5, aspect=5)
            plt.show()

        if "Flow" in names_of_fields:
            # FLOW
            fig1 = plt.figure(figsize=(8, 6), facecolor='w', edgecolor='k',
                              linewidth=2.0, frameon=True)
            ax = fig1.gca(projection='3d')
            surf = ax.plot_surface(self._T, self._X, self._q, rstride=4, cstride=4, alpha=0.2,
                                   cmap=cm.coolwarm,
                                   linewidth=0.4, antialiased=True)
            # ax.contourf(T, X, q * qc, zdir='y', offset=50, cmap=cm.viridis)
            # ax.plot_wireframe(T, X, q*qc, rstride=4, cstride=4, alpha=1.)
            ax.set_xlabel('t')
            ax.set_ylabel('z')
            ax.set_zlabel('flow')
            ax.set_title("Pulse flow wave")
            # ax.set_xlim([min(t_scaled[len(t_scaled)-(len(t_scaled)/4)::]), max(t_scaled[len(t_scaled)-(len(t_scaled)/4)::])])
            # ax.set_ylim([min(x_scaled), 60])
            cb = fig1.colorbar(surf, shrink=0.5, aspect=5)
            # cb.set_label("$mm^3/sec$")
            # plt.legend(loc='upper center', scatterpoints=1, ncol=3,
            #            bbox_to_anchor=(0.5, -0.03), fancybox=True,
            #            shadow=True)
            plt.show()

        if "Pressure" in names_of_fields:
            # PRESSURE
            fig2 = plt.figure(figsize=(8, 6), facecolor='w', edgecolor='k',
                              linewidth=2.0, frameon=True)
            ax = fig2.gca(projection='3d')
            # *7500.6156130264
            surf = ax.plot_surface(self._T, self._X, self._p,
                                   rstride=4, cstride=4, alpha=0.5,
                                   cmap=cm.coolwarm,
                                   linewidth=0.4, antialiased=True)
            # cset = ax.contourf(self._T, self._X, self._p,
            #                    zdir='y', offset=60, cmap=cm.viridis)
            # cset = ax.contourf(T, X, P, zdir='z', offset=300, cmap=cm.coolwarm, alpha=0.4)
            ax.set_xlabel('t')
            ax.set_ylabel('z')
            ax.set_zlabel('pressure')
            # ax.set_title(r"$p = p_{ext} + f(R_0, k)(1 - \sqrt{\frac{A_0}{A})}$")
            # ax.set_xlim([min(t_scaled[len(t_scaled)-(len(t_scaled)/4)::]), max(t_scaled[len(t_scaled)-(len(t_scaled)/4)::])])
            # ax.set_ylim([min(x_scaled), 60])
            # ax.set_zlim([0, 1700.])
            cb2 = fig2.colorbar(surf, shrink=0.5, aspect=5)
            # cb2.set_label("$MPa$")
            plt.show()

        if "Velocity" in names_of_fields:
            # VELOCITY
            fig4 = plt.figure(figsize=(8, 6), facecolor='w', edgecolor='k',
                              linewidth=2.0, frameon=True)
            ax = fig4.gca(projection='3d')
            surf = ax.plot_surface(self._T, self._X, self._u, rstride=4, cstride=4,
                                   alpha=0.25,
                                   linewidth=0.25, antialiased=True)
            # ax.view_init(30, -60)
            cset = ax.contourf(self._T, self._X, self._u, zdir='Z', offset=0,
                               cmap=cm.coolwarm, alpha=0.4, linewidth=0.3)
            csetl = ax.contour(self._T, self._X, self._u, zdir='Z', offset=5, linewidth=0.3)
            ax.set_xlabel('t')
            ax.set_ylabel('z')
            ax.set_zlabel('Velocity')
            ax.set_title(r"$\bar{u}$")
            # ax.set_xlim([min(t_scaled[len(t_scaled)-(len(t_scaled)/4)::])-1, max(t_scaled[len(t_scaled)-(len(t_scaled)/4)::])])
            # ax.set_ylim([min(x_scaled), max(x_scaled)])
            cb4 = fig4.colorbar(cset, shrink=0.5, aspect=5)
            # cb4.set_label("mm/sec")
            plt.show()

        if "PU-LOOP" in names_of_fields:
            # PU-LOOP
            fig6 = plt.figure(figsize=(8, 6), facecolor='w', edgecolor='k',
                              linewidth=2.0, frameon=True)
            ax = fig6.gca(projection='3d')
            surf = ax.plot_surface(self._X, self._u,
                                   self._p, rstride=4,
                                   cstride=4, alpha=0.25,
                                   linewidth=0.2, antialiased=True)
            ax.view_init(30, -40)
            # ax.plot(T[0, :], U[0, :], P[0, :], "k--", linewidth=2.)
            # for i in range(X.shape[0]):
            #     ax.plot(X[i, 0:5], U[i, 0:5], P[i, 0:5], "k-", linewidth=1.,
            #             label="wave speed")
            # ax.plot(T[-1, 0:8], U[-1, 0:8], P[-1, 0:8], "m-", linewidth=2.)

            cset = ax.contour(self._X, self._u,
                              self._p,
                              zdir='x', offset=0)
            ax.set_xlabel('x')
            ax.set_ylabel('velocity')
            ax.set_zlabel('pressure')
            ax.set_title('$PU$ Loop')
            # ax.set_xlim([min(x_scaled)-5, max(x_scaled)])
            # plt.legend(loc="upper right")
            # ax.set_ylim([min(x_scaled), max(x_scaled)])
            # fig6.colorbar(surf, shrink=0.5, aspect=5)
            plt.show()

    def plot_2Dfield(self, names_of_fields):
        if "QP-LOOP" in names_of_fields:
            # QP-LOOP
            rows, columns = self._p.shape
            fig_PU = plt.figure(figsize=(8, 6), facecolor='w', edgecolor='k',
                                linewidth=2.0, frameon=True)
            n_curves = 5
            step = (rows + 1) / n_curves
            for i in range(1, rows, step):
                plt.plot(self._q[i, :],
                         self._p[i, :],
                         label="QU loop point %d" % i)
            plt.xlabel("q")
            plt.ylabel("p")
            plt.title("QP-Loop")
            plt.legend(loc='upper center', scatterpoints=1, ncol=2,
                       bbox_to_anchor=(0.5, -0.03), fancybox=True,
                       shadow=True)
            plt.show()

    @property
    def T(self):
        return self._T

    @property
    def X(self):
        return self._X

    @property
    def u(self):
        return self._u

    @property
    def A(self):
        return self._A

    @property
    def p(self):
        return self._p

    @property
    def q(self):
        return self._q

if __name__ == '__main__':
    raise NotImplementedError('Module is not idented for direct execution')