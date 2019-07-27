"""
Pylsewave.mesh is a module that contains classes for vessel/tube and vessel network definition.
"""
from __future__ import division
import numpy as np

__author__ = "Georgios E. Ragkousis"
STATUS_OK = 0
STATUS_ERROR = -1


class Vessel(object):
    """
    Base class for singe-vessel/tube definition

    """

    def __init__(self, name, L, R_proximal, R_distal, Wall_thickness,
                 WindKessel=None, Id=0, model_type='Linear_elastic'):
        """
        Vessel-Class constructor

        :param name: name of the vessel/segment
        :type name: str
        :param L: Length of the segment
        :type L: float
        :param R_proximal: Proximal radius
        :type R_proximal: float
        :param R_distal: Distal radius
        :type R_distal: float
        :param Wall_thickness: (float) Wall thickness
        :type Wall_thickness: float
        :param WindKessel: Windkessel params (R, L, C)
        :type WindKessel: dict
        :type mapping: dict(str, float)
        :param Id: Unique id of the segment
        :type Id: int
        :param model_type: model of the vessel (linear-elastic, visco-elastic)
        :type model_type: str
        :param priority: "Linear Elastic" or "Visco-elastic"

        :Example:

        >>> myVessel = Vessel(name="Barchial artery", L=100, R_proximal=4.0,
                              R_disWindKessel=None, Id=0, model_type='Linear_elastic')

        """
        self._name = name
        self._L = L
        self._R_pro = R_proximal
        self._R_dist = R_distal
        self._R0 = None
        self._Wall_th = Wall_thickness
        self._Id = Id
        self._dx = None
        self._x = None
        self._RLC = None
        self._f_r0 = None
        self._df_dr0 = None
        self._df_dx = None
        self._f_r0_ph = None
        self._df_dr0_ph = None
        self._df_dx_ph = None
        self._f_r0_mh = None
        self._df_dr0_mh = None
        self._df_dx_mh = None
        self._k = None
        self._phi = None
        self._model_type = model_type
        # self._u = None
        # self._u_n = None
        if WindKessel is not None:
            if isinstance(WindKessel, dict):
                self._RLC = WindKessel
            else:
                raise ValueError("The Windkessel param should be type dict!")

    @property
    def id(self):
        """
        :getter: returns the unique vessel/segment Id
        :setter: sets the unique vessel/segment Id
        :type: (int)
        """
        return self._Id

    @id.setter
    def id(self, value):
        self._Id = value

    @property
    def name(self):
        """
        :getter: returns the name of the vessel
        :type: (string)
        """
        return self._name

    @property
    def length(self):
        """
        :getter: returns the Length of the 1D segment
        :type: float
        """
        return self._L

    @property
    def r_prox(self):
        """
        :getter: returns the: :math:`R_{proximal}`
        :type: float
        """
        return self._R_pro

    @property
    def r_dist(self):
        """
        :getter: returns the: :math:`R_{distal}`
        :type: float
        """
        return self._R_dist

    @property
    def w_th(self):
        """
        :getter: returns the :math:`W_{thickness}`
        :type: float
        """
        return self._Wall_th

    @property
    def dx(self):
        """
        The spatial descretisation of the Vessel object.

        .. math: `dx = \\frac{L}{Nx}`

        When we set this property/attribute of the object, several other Vessel
        characteristics are calculated, such as :py:meth:`mesh.Vessel.r0`,
        :py:meth:`mesh.Vessel.f_r0`, etc.

        :getter: returns the spatial dx of the vessel
        :setter: sets the dx of vessel/segment
        :type: (float)
        """
        return self._dx

    @property
    def phi(self):
        return self._phi

    @phi.setter
    def phi(self, value):
        self._phi = value

    @dx.setter
    def dx(self, dx):
        # here we reduce dx if the length of the segment is less than dx
        if (int(round(self.length/dx)) + 1) == 1:
            self._x = np.array([0., self.length])
        # self._dx = dx
        else:
            self._x = np.linspace(0., self.length, int(round(self.length/dx)) + 1)
        # make sure that dx is correct!
        self._dx = self._x[1] - self._x[0]
        # print self._dx
        # self._R0 = self.r_prox * np.exp(np.log(self.r_dist / self.r_prox) * (self.x / self.length))
        self._R0 = self.calculate_R0(self._x)
        if self._k is not None:
            self._f_r0 = self.f(self._R0, self._k)
            self._df_dr0 = self.dfdr(self._R0, self._k)
            self._df_dx = np.gradient(self._R0, self._dx)
            self._f_r0_ph  = self.f(self.interpolate_R0(0.5), self._k)
            self._df_dr0_ph = self.dfdr(self.interpolate_R0(0.5), self._k)
            self._df_dx_ph = np.gradient(self.interpolate_R0(0.5), self._dx)
            self._f_r0_mh  = self.f(self.interpolate_R0(-0.5), self._k)
            self._df_dr0_mh = self.dfdr(self.interpolate_R0(-0.5), self._k)
            self._df_dx_mh = np.gradient(self.interpolate_R0(-0.5), self._dx)

    @property
    def x(self):
        return self._x

    @property
    def RLC(self):
        """
        :getter: returns the RLC params
        :setter: sets the RLC params
        :type: dict
        :type mapping: dict(str, float)
        :raises ValueError: if the input is not a dictionary
        """
        return self._RLC

    @RLC.setter
    def RLC(self, rlc):
        if isinstance(rlc, dict):
            self._RLC = rlc
        else:
            raise ValueError("Input should be dictionary")

    @property
    def r0(self):
        """
        :getter: returns the reference radius along the 1D segment
        :type: ndarray[ndim=1, type=float]
        """
        return self._R0

    @property
    def f_r0(self):
        """
        This property contains the f(r0, k) function

        .. math: `f(R_0, k) = \\frac{4}{3}(k_1 e^{R_0 k_2} + k_3)`

        :getter: returns the f value calculated in every node of the 1D vessel
        """
        return self._f_r0

    @property
    def df_dr0(self):
        return self._df_dr0

    @property
    def df_dx(self):
        return self._df_dx

    @property
    def f_r0_ph(self):
        return self._f_r0_ph

    @property
    def df_dr0_ph(self):
        return self._df_dr0_ph

    @property
    def df_dx_ph(self):
        return self._df_dx_ph

    @property
    def f_r0_mh(self):
        return self._f_r0_mh

    @property
    def df_dr0_mh(self):
        return self._df_dr0_mh

    @property
    def df_dx_mh(self):
        return self._df_dx_mh

    @property
    def k(self):
        return self._k

    #added
    @staticmethod
    def beta(R0, E, h):
        return (4/3.)*(np.sqrt(np.pi)*E*h / (np.pi*R0*R0))

    @staticmethod
    def f(R0, k):
        """
        This property calculates the f(r0, k) function

        .. math: `f(R_0, k) = \\frac{4}{3}(k_1 e^{R_0 k_2} + k_3)`

        :param R0: Reference diameter
        :type priority: float or ndarray[ndim=1, type=float]
        :param k: the empirical k params
        :type k: ndarray[ndim=1, type=float]
        :return: float or ndarray[ndim=1, type=float]
        """
        k1, k2, k3 = k
        res = (4 / 3.) * (k2 * np.exp(k3 * R0) + k1)

        return res

    @staticmethod
    def dfdr(R0, k):
        """
        This property calculates the df(r0, k)/dr function

        .. math: `\\frac{f(R_0, k)}{dR_0} = \\frac{4}{3}(k_2 k_3 e^{R_0 k_3})`

        :param R0: Reference diameter
        :type priority: float or ndarray[ndim=1, type=float]
        :param k: the empirical k params
        :type k: ndarray[ndim=1, type=float]
        :return: float or ndarray[ndim=1, type=float]
        """
        k1, k2, k3 = k
        return (4 / 3.) * k2 * k3 * np.exp(k3 * R0)

    def set_k_vector(self, k):
        """
        Sets the k vector (k1, k2, k3) params
        :param k: the k vector
        :type k: ndarray[ndim=1, type=float]
        :return: None
        :raise ValueError: if the k input variable is not ndarray
        """
        if isinstance(k, np.ndarray):
            self._k = k
        else:
            raise ValueError("k vector should be type ndarray")

    def calculate_R0(self, x):
        """
        Method to calculate the reference radius :math:`R_0(x)`.
        The default model is

        .. math::
            R_0(x) = R_{prox} \\exp(\\log(R_{distal} / R_{prox})(x/L))

        :param x: (array or float) spatial point to calculate the reference diameter
        :return: (array or float) :math:`R_0(x)`
        """
        return self.r_prox*np.exp(np.log(self.r_dist/self.r_prox)*(x/self.length))

    def interpolate_R0(self, value):
        """
        Method to interpolate initial radius values in points offset from nodes

        :param value: the node offset interpolation points (e.g. 0.5, -0.5)
        :type value: float
        :return: float or ndarray[ndim=1, type=float]
        """
        return self.r_prox*np.exp(np.log(self.r_dist/self.r_prox)*((self.x + self.dx*value)/self.length))


class VesselNetwork(object):
    """
    Class to assemble the single vessels/segments that comprise the Vessel network
    """
    def __init__(self, vessels, dx=None, Nx=None):
        """
        :param vessels: (list) object containing the Vessel segments
        :param dx: (float) global spatial descretisation
        :param Nx: (int) global number of nodes per element
        :raises ValueError: if the vessels are not pylsewave.mesh.Vessel type
        """
        self.vessels = vessels
        self._bif = None
        self._conj = None

        if dx is not None:
            self._dx = dx
        if Nx is not None:
            self._Nx = Nx

        if isinstance(self.vessels, list):
            for i in range(len(self.vessels)):
                if isinstance(self.vessels[i], (Vessel)):
                    self.vessels[i].dx = self._dx
                else:
                    raise ValueError("The vessel no_%d (input) is not of"
                                     "Class Vessel/VesselScaled/cyVessel/cyVesselScaled!" % (i + 1))
        if Nx is not None:
            l_min = min([self.vessels[i].length for i in range(len(self.vessels))])
            Nx_min = int(round(l_min / dx))
            if Nx_min < Nx:
                self._dx = l_min / Nx
                for i in range(len(self.vessels)):
                    self.vessels[i].dx = self._dx

    @property
    def bif(self):
        """
        :getter: returns the vessel connectivity (bifurcations)
        :setter: sets the vessel connectivity (bifurcations)
        :type: ndarray[ndim=2, type=float]
        """
        return self._bif

    @bif.setter
    def bif(self, value):
        self._bif = value

    @property
    def conj(self):
        """
        :getter: returns the vessel connectivity (conjunctions)
        :setter: sets the vessel connectivity (conjunctions)
        :type: ndarray[ndim=2, type=float]
        """
        return self._conj

    @conj.setter
    def conj(self, value):
        self._conj = value

    @property
    def Nx(self):
        """
        :getter: returns the global Nx (number of nodes)
        :type: int
        """
        return self._Nx

    @property
    def dx(self):
        """
        :getter: returns the global dx
        :type: float
        """
        return self._dx

if __name__ == '__main__':
    raise NotImplementedError('Module is not idented for direct execution')