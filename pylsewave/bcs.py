"""
pylsewave module containing classes for boundary conditions.
"""
from __future__ import division

import numpy as np
from pylsewave.pdes import PDEm

__author__ = "Georgios E. Ragkousis"
__acknowledgement__ = "U_L (outlet bc function) implemented/recycled from vamPy (https://github.com/akdiem/vampy)"

class BCs(object):
    """
    Base class for BCs definition.
    """
    def __init__(self, inpdes, inletfun=None):
        """

        :param pylsewave.pdes.PDEm inpdes: the PDE class of the problem
        :param pylsewave.interpolate.CubicSpline.eval_spline inletfun: interpolation fun to calculate the inlet prescribed field
        :raises TypeError: if the inpdes is not pylsewave.pdes.PDEm type
        """
        if isinstance(inpdes, PDEm) is not True:
            raise TypeError("inpdes should be pylsewave.pdes.PDEm type")

        self.mesh = inpdes.mesh
        if inpdes._rho is not None:
            self._rho = inpdes._rho
        self._pdes = inpdes
        self.inlet_fun = inletfun
        self.vessels = inpdes.vessels

    @property
    def pdes(self):
        """
        :getter: returns the pde system
        :type: pylsewave.pdes.PDEm
        """
        return self._pdes

    def I(self, x, vessel_index=None):
        """
        Sets the initial condition of the state variables

        :param x: the spatial coordinate of the node
        :type priority: float or ndarray[ndim=1, type=float]
        :param int vessel_index: the unique Id of the vessel
        :return: ndarray[ndim=1, type=float]
        """
        if isinstance(x, np.ndarray) and (vessel_index is not None):
            ret = np.zeros((2, len(x)))
            return np.array([np.pi * self.vessels[vessel_index].r0 ** 2, 0.0 * x])
        else:
            return np.array([np.pi * self.vessels.r0 ** 2, 0.0 * x])

    def U_0(self, u, t, dx, dt, vessel_index, out):
        """
        Method to compute the inlet BCs.

        :Example:

          Let assume that we have a MacCormack scheme. Moreover, we solve the
          PDE system with respect to A, Q. We prescribe also flow waveform at
          the inlet :math:`Q_{x=0}(t)`. Then from continuity equation,
          the cross sectional area :math:`A_{x=0}^{t=n+1}` can be calculated as

          .. math::

              A_{x=0}^{t=n+1} = A_{x=0}^{t=n} - \\left(2 \\frac{dt}{dx} \\left( Q_{x=1}^{t=n+1} - Q_{x=0}^{t=n+1}\\right) \\right)

        :param u: the solution vector u_n
        :param t: the time step of the simulation
        :param dx: the spatial discretisation of the vessel
        :param dt: the space discretisation of the simulation
        :param vessel_index: the unique vessel index
        :param out: the solution array u
        :return: None
        """
        theta = dt / dx

        q_pres = self.inlet_fun(t)
        q_p1 = out[1, 1]


        A = out[0, 0] - (2. * theta) * (q_p1 - q_pres)

        out[0, 0] = A
        out[1, 0] = q_pres

    def U_L(self, u, t, dx, dt, vessel_index, out):
        """

        Class-method to compute the outlet BCs. This method implements the three-element
        Windkessel model (RCR) along with the iterative method proposed in
        Kolachalama et al 2007, *Predictive Haemodynamics in a One-Dimensional
        Carotid Artery Bifurcation. Part I Application to Stent Design*. In: IEEE Transactions
        on Biomedical Engineering. This method has also been implemented in
        `VamPy <https://github.com/akdiem/vampy>`_

        :param u: Solution vector u_n
        :param t: time step
        :param dx: vessel/segment dx
        :param dt: time discretisation dt
        :param vessel_index: unique vessel index
        :param out: solution vector u
        :return: None
        """
        _A, _q = u
        a_out = None
        q_out = None
        theta = dt / dx
        dt2 = dt / 2.
        p_out = self._pdes.pressure(_A[-1], -1, vessel_index)
        p0 = self._pdes.pressure(_A[-1], -1, vessel_index)

        q_m1 = out[1, -2]

        # a_out = A_n - theta * (q_out - u_m1[1]) # alternative
        k = 0
        R1 = self.vessels[vessel_index].RLC["R_1"]
        Rt = self.vessels[vessel_index].RLC["R_t"]
        Ct = self.vessels[vessel_index].RLC["C_t"]

        x = (dt / (R1 * Rt * Ct))
        while k < 1000:
            p_old = p0
            q_out = (x * p_out - x * (R1 + Rt) * _q[-1] +
                     (p0 - p_out)/R1 + _q[-1])

            a_out = _A[-1] - theta * (q_out - q_m1)
            p0 = self._pdes.pressure(a_out, index=-1, vessel_index=vessel_index)
            if abs(p_old - p0) < 1e-7:
                break
            k += 1
#        print k
        out[0, -1] = a_out
        out[1, -1] = q_out

    def bifurcation_R(self, x, u, dt, vessel_index_list):
        """
        Method to calculate the non-linear system of equations in the bifurcation.
        Usually, 6 equations are given.

        :param x: the solution vector x at an iteration i
        :param u: the solution vector u_n of internal nodes
        :param dt: time step dt
        :param vessel_index_list: unique vessel Id
        :return: ndarray[ndim=1, type=float]
        :raises NotImpementedError: this is abstract class, define an inherited class to override this method.
        """
        raise NotImplementedError

    def bifurcation_Jr(self, x, u, dt, vessel_index_list):
        """
        Method to calculate the Jacobian of the non-linear system formed at the bifurcations.
        Normally the output will be a 6x6 array/matrix.

        :param x: the solution vector at iteration i
        :param u: the solution vector u_n of internal nodes
        :param dt: time step dt
        :param vessel_index_list: the unique vessel Id
        :return: ndarray[ndim=2, dtype=float]
        :raises NotImpementedError: this is abstract class, define an inherited class to override this method.
        """
        raise NotImplementedError

    def conjuction_R(self, x, u, dt, vessel_index_list):
        """
        Method to calculate the non-linear system of equations in a conjunction.
        Usually, 4 equations are given.

        :param x: the solution vector x at an iteration i
        :param u: the solution vector u_n of internal nodes
        :param dt: time step dt
        :param vessel_index_list: unique vessel Id
        :return: ndarray[ndim=1, type=float]
        :raises NotImpementedError: this is abstract class, define an inherited class to override this method.
        """
        raise NotImplementedError

    def conjuction_Jr(self, x, u, dt, vessel_index_list):
        """
        Method to calculate the Jacobian of the non-linear system formed at the conjunctions.
        Normally the output will be a 4x4 array/matrix.

        :param x: the solution vector at iteration i
        :param u: the solution vector u_n of internal nodes
        :param dt: time step dt
        :param vessel_index_list: the unique vessel Id
        :return: ndarray[ndim=2, dtype=float]
        :raises NotImpementedError: this is abstract class, define an inherited class to override this method.
        """
        raise NotImplementedError


class BCsWat(BCs):
    """
    Class for BCs definition.
    pylsewave.bcs.BCs inherited class for the pylsewave.pdes.PDEsWat class.
    """
    def U_0(self, u, t, dx, dt, vessel_index, out):
        """
        Method to compute the inlet BCs.
        This method assumes that a pressure is prescribed at the inlet.

        .. note::

            override the method in case you want to prescribe other type of inlet BCs.

        :param u: the solution vector u_n
        :param t: the time step of the simulation
        :param dx: the spatial discretisation of the vessel
        :param dt: the space discretisation of the simulation
        :param vessel_index: the unique vessel index
        :param out: the solution array u
        :return: None
        """
        _A, _q = u
        theta = dt / dx
        dt2 = dt / 2.
        p_pres = self.inlet_fun(t)
        if vessel_index is None:
            A0 = np.pi*self.vessels.r0[0]**2.0
            f = self.vessels.f_r0[0]
        else:
            A0 = np.pi*self.vessels[vessel_index].r0[0]**2.0
            f = self.vessels[vessel_index].f_r0[0]

        A = A0*((p_pres /f) + 1)**2.0
        W2_0 = _q[0]/_A[0] - 4*self._pdes.compute_c(_A[0], index=0, vessel_index=vessel_index)
        W2_1 = _q[1]/_A[1] - 4*self._pdes.compute_c(_A[1], index=1, vessel_index=vessel_index)
        l2 = _q[0]/_A[0] - self._pdes.compute_c(_A[0], index=0, vessel_index=vessel_index)
        x = -l2*dt
#         print x
#         print("x", x)
        # W2 = _q[1]/_A[1] - 4*self.pdes.compute_c_i(_A[1], index=1, vessel_index=vessel_index)
        W2 = self.linear_extrapolation(x, 0.0, dx, W2_0, W2_1)
#         print("W2_0", W2_0)
#         print("W2_1", W2_1)
#         print("W2", W2)

        q = A*(W2 + 4*self._pdes.compute_c(A, index=0, vessel_index=vessel_index))
#         print A
#         print q
        out[0, 0] = A
        out[1, 0] = q

    @staticmethod
    def linear_extrapolation(x, x1, x2, y1, y2):
        """
        Static method for linear extrapolation..

        :param x: location of extrapolation
        :param x1: x value left
        :param x2: y value left
        :param y1: x value right
        :param y2: y value right
        :return: extrapolation value
        """
        return y1 + (x - x1)*((y2 - y1)/(x2-x1))

    def bifurcation_R(self, x, u, dt, vessel_index_list):
        """
        Method to calculate the non-linear system of equations in the bifurcation.
        Usually, 6 equations are given.

        :param x: the solution vector x at an iteration i
        :param u: the solution vector u_n of internal nodes
        :param dt: time step dt
        :param vessel_index_list: unique vessel Id
        :return: ndarray[ndim=1, type=float]
        """
        x1, x2, x3, x4, x5, x6 = x
        p, d1, d2 = vessel_index_list
        A_p, q_p, A_d1, q_d1, A_d2, q_d2 = u

        return np.array([-x2 + x4 + x6,
                         (-0.5*self._rho*(x2/x1)**2 - self.vessels[p].f_r0[-1]*(
                         np.sqrt(x1 / (np.pi*self.vessels[p].r0[-1]**2)) - 1) + 0.5*self._rho*(x4/x3)**2 +
                          self.vessels[d1].f_r0[0]*(np.sqrt(x3 / (np.pi*self.vessels[d1].r0[0]**2)) - 1)),
                         (-0.5*self._rho*(x2/x1)**2 - self.vessels[p].f_r0[-1]*(
                         np.sqrt(x1 / (np.pi*self.vessels[p].r0[-1]**2)) - 1) + 0.5*self._rho*(x6/x5)**2 +
                          self.vessels[d2].f_r0[0]*(np.sqrt(x5 / (np.pi*self.vessels[d2].r0[0]**2)) - 1)),
                         (-x2/x1 - 4*self._pdes.compute_c(x1, index=-1, vessel_index=p) +
                          (q_p/A_p + 4*self._pdes.compute_c(A_p, index=-2, vessel_index=p))),
                         (-x4/x3 + 4*self._pdes.compute_c(x3, index=0, vessel_index=d1) +
                          (q_d1/A_d1 - 4*self._pdes.compute_c(A_d1, index=1, vessel_index=d1))),
                         (-x6/x5 + 4*self._pdes.compute_c(x5, index=0, vessel_index=d2) +
                          (q_d2/A_d2 - 4*self._pdes.compute_c(A_d2, index=1, vessel_index=d2)))])

    def bifurcation_Jr(self, x, u, dt, vessel_index_list):
        """
         Method to calculate the Jacobian of the non-linear system formed at the bifurcations.
         Normally the output will be a 6x6 array/matrix.

         :param x: the solution vector at iteration i
         :param u: the solution vector u_n of internal nodes
         :param dt: time step dt
         :param vessel_index_list: the unique vessel Id
         :return: ndarray[ndim=2, dtype=float]
         """
        x1, x2, x3, x4, x5, x6 = x
        p, d1, d2 = vessel_index_list
        return np.array([[0., -1., 0., 1., 0., 1.],
                         [self._rho*(x2**2 / x1**3) - 0.5*self.vessels[p].f_r0[-1]*((np.pi*self.vessels[p].r0[-1]**2)**(-0.5))*(x1**(-1/2.)),
                          -self._rho*x2 / x1**2,
                          -self._rho*(x4**2 / x3**3) + 0.5*self.vessels[d1].f_r0[0]*((np.pi*self.vessels[d1].r0[0]**2)**(-0.5))*(x3**(-1/2.)),
                          self._rho*x4 / x3**2, 0., 0.],
                         [self._rho*(x2**2 / x1**3) - 0.5*self.vessels[p].f_r0[-1]*((np.pi*self.vessels[p].r0[-1]**2)**(-0.5))*(x1**(-1/2.)),
                          -self._rho*x2 / x1**2, 0., 0.,
                          -self._rho*(x6**2 / x5**3) + 0.5*self.vessels[d2].f_r0[0]*((np.pi*self.vessels[d2].r0[0]**2)**(-0.5))*(x5**(-1/2.)),
                          self._rho*x6 / x5**2],
                         [(x2 / x1**2) - 4*(1/4.*np.sqrt((0.5/self._rho)*self.vessels[p].f_r0[-1] / np.sqrt(np.pi*self.vessels[p].r0[-1]**2)))*x1**(-3/4.),
                          -1./x1, 0., 0., 0., 0.],
                         [0., 0.,
                          (x4 / x3**2) + 4*(1/4.*np.sqrt((0.5/self._rho)*self.vessels[d1].f_r0[0] / np.sqrt(np.pi*self.vessels[d1].r0[0]**2)))*x3**(-3/4.),
                          -1./x3, 0., 0.],
                         [0., 0., 0., 0.,
                          (x6 / x5**2) + 4*(1/4.*np.sqrt((0.5/self._rho)*self.vessels[d2].f_r0[0] / np.sqrt(np.pi*self.vessels[d2].r0[0]**2)))*x5**(-3/4.),
                          -1./x5]])

    def conjuction_R(self, x, u, dt, vessel_index_list):
        """
        Method to calculate the non-linear system of equations in a conjunction.
        Usually, 4 equations are given.

        :param x: the solution vector x at an iteration i
        :param u: the solution vector u_n of internal nodes
        :param dt: time step dt
        :param vessel_index_list: unique vessel Id
        :return: ndarray[ndim=1, type=float]
        """
        x1, x2, x3, x4 = x
        d1, d2 = vessel_index_list
        A_d1, q_d1, A_d2, q_d2 = u

        return np.array([-x2 + x4,
                         (-0.5*self._rho*(x2/x1)**2 - self.vessels[d1].f_r0[-1]*(
                         np.sqrt(x1 / (np.pi*self.vessels[d1].r0[-1]**2)) - 1) + 0.5*self._rho*(x4/x3)**2 +
                         self.vessels[d2].f_r0[0]*(np.sqrt(x3/ (np.pi*self.vessels[d2].r0[0]**2)) - 1)),
                         (-x2/x1 - 4*self._pdes.compute_c(x1, index=-1, vessel_index=d1) +
                         (q_d1/A_d1 + 4*self._pdes.compute_c(A_d1, index=-2, vessel_index=d1))),
                         (-x4/x3 + 4*self._pdes.compute_c(x3, index=0, vessel_index=d2) +
                         (q_d2/A_d2 - 4*self._pdes.compute_c(A_d2, index=1, vessel_index=d2)))])

    def conjuction_Jr(self, x, u, dt, vessel_index_list):
        """
        Method to calculate the Jacobian of the non-linear system formed at the conjunctions.
        Normally the output will be a 4x4 array/matrix.

        :param x: the solution vector at iteration i
        :param u: the solution vector u_n of internal nodes
        :param dt: time step dt
        :param vessel_index_list: the unique vessel Id
        :return: ndarray[ndim=2, dtype=float]
        """
        x1, x2, x3, x4 = x
        d1, d2 = vessel_index_list

        return np.array([[0., -1., 0., 1.],
                         [self._rho*(x2**2 / x1**3) - 0.5*self.vessels[d1].f_r0[-1]*((np.pi*self.vessels[d1].r0[-1]**2)**(-0.5))*(x1**(-1/2.)),
                          -self._rho*x2 / x1**2,
                          -self._rho*(x4**2 / x3**3) + 0.5*self.vessels[d2].f_r0[0]*((np.pi*self.vessels[d2].r0[0]**2)**(-0.5))*(x3**(-1/2.)),
                          self._rho*x4 / x3**2],
                         [(x2 / x1**2) - 4*(1/4.*np.sqrt((0.5/self._rho)*self.vessels[d1].f_r0[-1] / np.sqrt(np.pi*self.vessels[d1].r0[-1]**2)))*x1**(-3/4.),
                          -1./x1, 0., 0.],
                         [0., 0.,
                          (x4 / x3**2) + 4*(1/4.*np.sqrt((0.5/self._rho)*self.vessels[d2].f_r0[0] / np.sqrt(np.pi*self.vessels[d2].r0[0]**2)))*x3**(-3/4.),
                          -1./x3]])

    if __name__ == '__main__':
        raise NotImplementedError('Module is not idented for direct execution')