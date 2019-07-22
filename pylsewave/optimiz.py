import numpy as np
from scipy.optimize import least_squares, differential_evolution

__author__ = "Georgios E. Ragkousis"

class OptimiseKVector(object):
    """
    Class for creating optimiser (ga or least-square)
    """
    def __init__(self, data, f, jac, method="least_squares"):
        """
        Class for optimising k elasticity parameters with least
         squares method or with a ga

        :param data: (float[]) data to be fitted by f function
        :param f: (void) function to be fitted to the given data
        :param jac: (void) function to compute the Jacobian of f
        :param method: (string) parameter to define the optimisation method
        :return: (float[]) k vector with optimised parameters
        """
        self._data = data
        self._f = f
        self._jac = jac
        self._method = method
        self._k = None
        self._k0 = None

    def _residual_f(self, x, t, y):
        return self.f(x, t) - y

    def _residual_f_ga(self, x, *data):
        t, y = data
        res = 0
        for i in range(len(t)):
            res += (self.f(x, t[i]) - y[i])**2.

        return res*0.5

    def _func(self, parameters, *data):
        # we have 3 parameters which will be passed as parameters and
        # "experimental" x,y which will be passed as data
        a, b, c = parameters
        x, y = data

        result = 0

        for i in range(len(x)):
            result += (a*x[i]**2 + b*x[i] + c - y[i])**2

        return result**0.5

    def set_initial_values(self, k0):
        if isinstance(k0, np.ndarray) is False or np.ndim(k0) > 1:
            raise ValueError("k0 should be 1D numpy array")
        else:
            self._k0 = k0

    def optimize_k(self, loss="cauchy", verbose=0):
        """
        method to optimise k values

        :param loss: type of loss function (see scipy.optimise.least_squares)
        :param verbose: print the output
        :return: None
        """
        _args = (self.data[:, 0], self.data[:, 1])
        _fit = least_squares(self._residual_f, self.k0, jac=self.jac, loss="cauchy",
                             f_scale=0.1, verbose=verbose, args=_args)
        self._k = _fit.x

    @property
    def data(self):
        return self._data

    @property
    def f(self):
        return self._f

    @property
    def jac(self):
        return self._jac

    @property
    def method(self):
        return  self._method

    @property
    def k0(self):
        return self._k0

    @property
    def k(self):
        return self._k