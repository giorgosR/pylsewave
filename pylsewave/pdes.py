__author__ = "Georgios E. Ragkousis"
import numpy as np
from pwutils import h_walls
from mesh import VesselNetwork

class PDEm(object):
    """
    Base class for PDE system declaration.
    """
    def __init__(self, mesh, rho, mu, Re=None):
        """
        .. math:

            \\frac{\\partial U}{\\partial t} + \\frac{\\partial F}{\\partial x} = S(U)

        :param pylsewave.mesh.VesselNetwork mesh: the vessel network
        :param float rho: the density of fluid medium (e.g. blood, water, etc.)
        :param float Re: the Reynolds number
        :raises TypeError: if the mesh is not pylsewave.mesh.VesselNetwork
        """
        if isinstance(mesh, VesselNetwork) is not True:
            raise TypeError("mesh should be of VesselNetwork type!")
        self.mesh = mesh
        if Re is not None:
            self._Re = Re
        if rho is not None:
            self._rho = rho
        if mu is not None:
            self._mu = mu
        self.vessels = mesh.vessels
        self._delta = None

    def set_boundary_layer_th(self, T, no_cycles):
        """
        method to calculate the boundary layer

        .. math::

            \\delta = \\sqrt{\\frac{\\nu T_{cycle}}{2 \\pi}}

        :param T: (float) Total period of the simulation
        :param no_cycles: No of cycles that the simulation will run
        :return: (float) Boundary layer delta
        """
        _nu = self._mu / self._rho
        T_cycle = float(T) / no_cycles
        self._delta = np.power((_nu * T_cycle) / (2. * np.pi), 0.5)

    @property
    def delta(self):
        return self._delta


    def compute_c(self, A, index=None, vessel_index=None):
        """
        Method to compute local wave speed

        :param A: cross-sectional area of vessel
        :type priority: float or ndarray[ndim=1, type=float]
        :param int index: the number or the node to be calculated
        :param vessel_index: the vessel unique Id
        :return: float or ndarray[ndim=1, type=float]
        """

        raise NotImplementedError

    def pressure(self, a, index=None, vessel_index=None):
        """
        Method to compute pressure from pressure-area relationship

        :param a: cross-sectional area of the vessel
        :type priority: float or ndarray[ndim=1, type=float]
        :param int index: index of the node in vessel
        :param int vessel_index: the vessel unique Id
        :return: float or ndarray[ndim=1, type=float]
        """
        raise NotImplementedError

    def flux(self, u, x, index=None, vessel_index=None):
        """
        Method to calculate the Flux term.

        :param u: the solution vector u
        :type priority: float or ndarray[ndim=2, dtype=float]
        :param x: the spatial points discretising the 1D segments
        :type priority: float or ndarray[ndim=2, dtype=float]
        :param int index: the index Id of the node
        :param int vessel_index: the unique vessel Id
        :return: ndarray[ndim=1, dtype=float] or ndarray[ndim=2, dtype=float]
        :raises NotImplementedError: This is a base class, define an inherited class to override this method
        """
        raise NotImplementedError

    def source(self, u, x, index=None, indexVessel=None):
        """
        Method to calculate the Source term.

        :param u: the solution vector u
        :type priority: float or ndarray[ndim=2, dtype=float]
        :param x: the spatial points discretising the 1D segments
        :type priority: float or ndarray[ndim=2, dtype=float]
        :param int index: the index Id of the node
        :param int vessel_index: the unique vessel Id
        :return: ndarray[ndim=1, dtype=float] or ndarray[ndim=2, dtype=float]
        :raises NotImplementedError: This is a base class, define an inherited class to override this method
        """
        raise NotImplementedError


class PDEsWat(PDEm):
    """
    Class for PDE system as defined in Watanabe et al. 2013, *Mathematical model of
    blood flow in an anatomically detailed arterial network of the arm*, In: Mathematical
    modelling and numerical analysis.
    """
    def compute_c(self, A, index=None, vessel_index=None):
        """
        Method to compute local wave speed

        .. math::
             c(A) = \\sqrt{\\frac{1}{2 \\rho} f(R_0, k) \\sqrt{\\frac{A}{A_0}}}

        :param A: cross-sectional area of vessel
        :type priority: float or ndarray[ndim=1, type=float]
        :param int index: the number or the node to be calculated
        :param vessel_index: the vessel unique Id
        :return: float or ndarray[ndim=1, type=float]
        """
        if vessel_index is not None:
            A0 = np.pi*self.vessels[vessel_index].r0**2
            if index is None:
                return np.sqrt((0.5/self._rho)*self.vessels[vessel_index].f_r0*np.sqrt(A / A0))
            else:
                return np.sqrt((0.5/self._rho)*self.vessels[vessel_index].f_r0[index]*np.sqrt(A / A0[index]))
        else:
            A0 = np.pi*self.vessels.r0**2
            if index is None:
                return np.sqrt((0.5/self._rho)*self.vessels.f_r0*np.sqrt(A / A0))
            else:
                return np.sqrt((0.5/self._rho)*self.vessels.f_r0[index]*np.sqrt( A / A0[index]))

    def pressure(self, a, index=None, vessel_index=None):
        """
        Method to compute pressure from pressure-area relationship

        .. math::
             p(A) = f(R_0, k) \\left( \\sqrt{\\frac{A}{A_0}} - 1 \\right)

        :param a: cross-sectional area of the vessel
        :type priority: float or ndarray[ndim=1, type=float]
        :param int index: index of the node in vessel
        :param int vessel_index: the vessel unique Id
        :return: float or ndarray[ndim=1, type=float]
        """
        if np.ndim(a) > 1 and isinstance(a, np.ndarray):
            raise ValueError("input array should be one dimensional")
        if vessel_index is not None:
            R0 = self.vessels[vessel_index].r0
            A0 = np.pi*R0**2
            if index == None:
                return self.vessels[vessel_index].f_r0*(np.sqrt(a / A0) - 1)
            else:
                return self.vessels[vessel_index].f_r0[index]*(np.sqrt(a / A0[index]) - 1)
        else:
            R0 = self.vessels.r0
            A0 = np.pi*R0**2
            if index is None:
                return self.vessels.f_r0*(np.sqrt(a / A0) - 1)
            else:
                return self.vessels.f_r0[index]*(np.sqrt(a / A0[index]) - 1)

    def flux(self, u, x, index=None, vessel_index=None):
        """
        Method to calculate the Flux term.

        :param u: the solution vector u
        :type priority: float or ndarray[ndim=2, dtype=float]
        :param x: the spatial points discretising the 1D segments
        :type priority: float or ndarray[ndim=2, dtype=float]
        :param int index: the index Id of the node
        :param int vessel_index: the unique vessel Id
        :return: ndarray[ndim=1, dtype=float] or ndarray[ndim=2, dtype=float]
        """
        _A, _q = u
        if np.ndim(u) > 1 and vessel_index is None:
            if x == 0:
                A0 = np.pi*self.vessels.r0**2
                f_r0 = self.vessels.f_r0
            elif x == 1:
                A0 = np.pi*(self.vessels.interpolate_R0(0.5))**2
                f_r0 = self.vessels.f_r0_ph
            elif x == -1:
                A0 = np.pi*(self.vessels.interpolate_R0(-0.5))**2
                f_r0 = self.vessels.f_r0_mh

            arr_out = np.zeros((2, len(self.vessels.x)))
            arr_out[0, :] = _q
            arr_out[1, :] = ((_q * _q) / _A) + (f_r0/(3.*self._rho))*((_A**(3/2.))*(A0**(-1/2.)))
            if index is None:
                return arr_out
            else:
                return arr_out[:, index]

        elif (np.ndim(u) > 1) and (vessel_index is not None):
            if x == 0:
                A0 = np.pi*self.vessels[vessel_index].r0**2
                f_r0 = self.vessels[vessel_index].f_r0
            elif x == 1:
                A0 = np.pi*(self.vessels[vessel_index].interpolate_R0(0.5))**2
                f_r0 = self.vessels[vessel_index].f_r0_ph
            elif x == -1:
                A0 = np.pi*(self.vessels[vessel_index].interpolate_R0(-0.5))**2
                f_r0 = self.vessels[vessel_index].f_r0_mh

            arr_out = np.zeros((2, len(self.vessels[vessel_index].x)))
            arr_out[0, :] = _q
            arr_out[1, :] = ((_q * _q) / _A) + (f_r0/(3*self._rho))*((_A**(3/2.))*(A0**(-1/2.)))
            if index is None:
                return arr_out
            else:
                return arr_out[:, index]

        elif (np.ndim(u) == 1) and (vessel_index is None):
            if x == 0:
                A0 = np.pi * self.vessels.r0 ** 2
                f_r0 = self.vessels.f_r0
            elif x == 1:
                A0 = np.pi * (self.vessels.interpolate_R0(0.5)) ** 2
                f_r0 = self.vessels.f_r0_ph
            elif x == -1:
                A0 = np.pi * (self.vessels.interpolate_R0(-0.5)) ** 2
                f_r0 = self.vessels.f_r0_mh
            if index is None:
                return np.array([_q, ((_q * _q) / _A) + (f_r0/(3*self._rho))*((_A**(3/2.))*(A0**(-1/2.)))])
            else:
                return np.array([_q, ((_q * _q) / _A) + (f_r0[index]/(3*self._rho))*((_A**(3/2.))*(A0[index]**(-1/2.)))])

        elif (np.ndim(u) == 1) and (vessel_index is not None):
            if x == 0:
                A0 = np.pi * self.vessels[vessel_index].r0 ** 2
                f_r0 = self.vessels[vessel_index].f_r0
            elif x == 1:
                A0 = np.pi * (self.vessels[vessel_index].interpolate_R0(0.5)) ** 2
                f_r0 = self.vessels[vessel_index].f_r0_ph
            elif x == -1:
                A0 = np.pi * (self.vessels[vessel_index].interpolate_R0(-0.5)) ** 2
                f_r0 = self.vessels[vessel_index].f_r0_mh
            if index is None:
                return np.array([_q, ((_q * _q) / _A) + (f_r0/(3*self._rho))*((_A**(3/2.))*(A0**(-1/2.)))])
            else:
                return np.array([_q, ((_q * _q) / _A) + (f_r0[index]/(3*self._rho))*((_A**(3/2.))*(A0[index]**(-1/2.)))])

    def source(self, u, x, index=None, indexVessel=None):
        """
        Method to calculate the Source term.

        :param u: the solution vector u
        :type priority: float or ndarray[ndim=2, dtype=float]
        :param x: the spatial points discretising the 1D segments
        :type priority: float or ndarray[ndim=2, dtype=float]
        :param int index: the index Id of the node
        :param int vessel_index: the unique vessel Id
        :return: ndarray[ndim=1, dtype=float] or ndarray[ndim=2, dtype=float]
        """
        _A, _q = u
        R = np.sqrt(_A / np.pi)
        _nu = (self._mu)/self._rho
        if np.ndim(u) > 1 and (indexVessel is None):
            if x == 0:
                A0 = np.pi * self.vessels.r0 ** 2
                f_r0 = self.vessels.f_r0
                df_dr0 = self.vessels.df_dr0
                df_dx = self.vessels.df_dx
            elif x == 1:
                A0 = np.pi * (self.vessels.interpolate_R0(0.5)) ** 2
                f_r0 = self.vessels.f_r0_ph
                df_dr0 = self.vessels.df_dr0_ph
                df_dx = self.vessels.df_dx_ph
            elif x == -1:
                A0 = np.pi * (self.vessels.interpolate_R0(-0.5)) ** 2
                f_r0 = self.vessels.f_r0_mh
                df_dr0 = self.vessels.df_dr0_mh
                df_dx = self.vessels.df_dx_mh

            arr_out = np.zeros((2, len(self.vessels.x)))
            A = (-2.*np.pi*_nu*R*_q) / (self._delta*_A)
            B = (1./self._rho)*(((2*np.sqrt(np.pi)*f_r0*_A**(3/2.))/(3.*A0)) - (((2/3.)*(_A**(3/2.))*A0**(-1/2.)) - _A)*df_dr0)
            C = df_dx
            arr_out[1, :] = A[:] + B[:] * C[:]
            if index is None:
                return arr_out
            else:
                return arr_out[:, index]

        elif (np.ndim(u) > 1) and (indexVessel is not None):
            if x == 0:
                A0 = np.pi * self.vessels[indexVessel].r0 ** 2
                f_r0 = self.vessels[indexVessel].f_r0
                df_dr0 = self.vessels[indexVessel].df_dr0
                df_dx = self.vessels[indexVessel].df_dx
            elif x == 1:
                A0 = np.pi * (self.vessels[indexVessel].interpolate_R0(0.5)) ** 2
                f_r0 = self.vessels[indexVessel].f_r0_ph
                df_dr0 = self.vessels[indexVessel].df_dr0_ph
                df_dx = self.vessels[indexVessel].df_dx_ph
            elif x == -1:
                A0 = np.pi * (self.vessels[indexVessel].interpolate_R0(-0.5)) ** 2
                f_r0 = self.vessels[indexVessel].f_r0_mh
                df_dr0 = self.vessels[indexVessel].df_dr0_mh
                df_dx = self.vessels[indexVessel].df_dx_mh

            arr_out = np.zeros((2, len(self.vessels[indexVessel].x)))
            A = (-2.*np.pi*_nu*R*_q) / (self._delta*_A)
            B = (1./self._rho)*(((2*np.sqrt(np.pi)*f_r0*_A**(3/2.))/(3.*A0)) - (((2/3.)*(_A**(3/2.))*A0**(-1/2.)) - _A)*df_dr0)
            C = df_dx
            arr_out[1, :] = A[:] + B[:] * C[:]
            if index is None:
                return arr_out
            else:
                return arr_out[:, index]

        elif (np.ndim(u) == 1) and (indexVessel is None):
            if x == 0:
                A0 = np.pi * self.vessels.r0 ** 2
                f_r0 = self.vessels.f_r0
                df_dr0 = self.vessels.df_dr0
                df_dx = self.vessels.df_dx
            elif x == 1:
                A0 = np.pi * (self.vessels.interpolate_R0(0.5)) ** 2
                f_r0 = self.vessels.f_r0_ph
                df_dr0 = self.vessels.df_dr0_ph
                df_dx = self.vessels.df_dx_ph
            elif x == -1:
                A0 = np.pi * (self.vessels.interpolate_R0(-0.5)) ** 2
                f_r0 = self.vessels.f_r0_mh
                df_dr0 = self.vessels.df_dr0_mh
                df_dx = self.vessels.df_dx_mh

            A = (-2.*np.pi*_nu*R*_q) / (self._delta*_A)
            B = (1./self._rho)*(((2*np.sqrt(np.pi)*f_r0*_A**(3/2.))/(3.*A0)) - (((2/3.)*(_A**(3/2.))*A0**(-1/2.)) - _A)*df_dr0)
            C = df_dx
            if index is None:
                return np.array([0.0, A + B * C])
            else:
                return np.array([0.0, A + B[index] * C[index]])

        elif (np.ndim(u) == 1) and (indexVessel is not None):
            if x == 0:
                A0 = np.pi * self.vessels[indexVessel].r0 ** 2
                f_r0 = self.vessels[indexVessel].f_r0
                df_dr0 = self.vessels[indexVessel].df_dr0
                df_dx = self.vessels[indexVessel].df_dx
            elif x == 1:
                A0 = np.pi * (self.vessels[indexVessel].interpolate_R0(0.5)) ** 2
                f_r0 = self.vessels[indexVessel].f_r0_ph
                df_dr0 = self.vessels[indexVessel].df_dr0_ph
                df_dx = self.vessels[indexVessel].df_dx_ph
            elif x == -1:
                A0 = np.pi * (self.vessels[indexVessel].interpolate_R0(-0.5)) ** 2
                f_r0 = self.vessels[indexVessel].f_r0_mh
                df_dr0 = self.vessels[indexVessel].df_dr0_mh
                df_dx = self.vessels[indexVessel].df_dx_mh

            A = (-2.*np.pi*_nu*R*_q) / (self._delta*_A)
            B = (1./self._rho)*(((2*np.sqrt(np.pi)*f_r0*_A**(3/2.))/(3.*A0)) - (((2/3.)*(_A**(3/2.))*A0**(-1/2.)) - _A)*df_dr0)
            C = df_dx
            if index is None:
                return np.array([0.0, A + B * C])
            else:
                return np.array([0.0, A + B[index] * C[index]])

class PDEsVisco(PDEsWat):
    """
    PDEs class inherited from PDEsWat. This PDE system contains the viscous part
    of the model, as well.
    """

    def cv_f_qua(self, a, w_th, vessel_index):
        A0 = np.pi * self.vessels[vessel_index].r0 ** 2
        return (2. / 3) * ((np.sqrt(np.pi) *self.vessels[vessel_index].phi*w_th) / (A0*np.sqrt(a)))

    def pressure(self, a, q, index=None, vessel_index=None):
        w_th = h_walls(self.vessels[vessel_index].r0)
        C_v = self.cv_f_qua(a, w_th, vessel_index)
        dqdx = np.gradient(q, self.vessels[vessel_index].dx)

        R0 = self.vessels[vessel_index].r0
        A0 = np.pi * R0 ** 2
        if index == None:
            return self.vessels[vessel_index].f_r0*(np.sqrt(a / A0) - 1) - C_v*dqdx
#            return self.vessels[vessel_index].f_r0 * (np.sqrt(a / A0) - 1) - cv * dqdx
        else:
            return self.vessels[vessel_index].f_r0[index]*(np.sqrt(a[index] / A0[index]) - 1) - C_v[index]*dqdx[index]
#            return self.vessels[vessel_index].f_r0[index] * (np.sqrt(a[index] / A0[index]) - 1) - cv[index] * dqdx[index]

    def pressure_i(self, a, q, index=None, vessel_index=None):
        R0 = self.vessels[vessel_index].r0
        A0 = np.pi * R0 ** 2
        if index == None:
            return self.vessels[vessel_index].f_r0*(np.sqrt(a / A0) - 1)
#            return self.vessels[vessel_index].f_r0 * (np.sqrt(a / A0) - 1) - cv * dqdx
        else:
            return self.vessels[vessel_index].f_r0[index]*(np.sqrt(a[index] / A0[index]) - 1)

if __name__ == '__main__':
    raise NotImplementedError('Module is not indented for direct execution')