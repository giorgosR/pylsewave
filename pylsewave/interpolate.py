import numpy as np

__author__ = "Georgios E. Ragkousis"


class CubicSpline(object):
    """
    Class for creating a cubic spline interpolator
    """
    def __init__(self, xdata, ydata):
        """
        Constructor

        :param xdata: array or list with x data
        :param ydata: array or list with y data
        """
        self._xdata = xdata
        self._ydata = ydata
        self._k = self._compute_curvatures()

    @staticmethod
    def ludecomp3(c, d, e):
        """
        c, d, e = LUdecomp3(c, d, e)
        LU decomposition of tridiagonal matrix [c\d\e]. On output, {c}, {d} and
        {e} are the diagonals of the decomposed matrix.
        """
        n = len(d)
        for k in range(1, n):
            lam = c[k - 1] / d[k - 1]
            d[k] = d[k] - lam * e[k - 1]
            c[k - 1] = lam
        return c, d, e

    @staticmethod
    def lusolve3(c, d, e, b):
        """
        x = LUsolve3(c, d, e, b)
        solves [c\d\e]{x} = {b}, where {c}, {d} and {e} are the vectors returned
        from LUdecomp3.
        """
        n = len(d)
        for k in range(1, n):
            b[k] = b[k] - c[k - 1] * b[k - 1]
        b[n - 1] = b[n - 1] / d[n - 1]
        for k in range(n - 2, -1, -1):
            b[k] = (b[k] - e[k] * b[k + 1]) / d[k]
        return b

    def _compute_curvatures(self):
        """
        k = curvatures(xData, yData)
        Returns the curvatures of cubic spline at its knots.
        """
        n = len(self.xdata) - 1
        c = np.zeros(n)
        d = np.ones(n + 1)
        e = np.zeros(n)
        k = np.zeros(n + 1)
        c[0:n - 1] = self.xdata[0:n - 1] - self.xdata[1:n]
        d[1:n] = 2.0 * (self.xdata[0:n - 1] - self.xdata[2:n + 1])
        e[1:n] = self.xdata[1:n] - self.xdata[2:n + 1]
        k[1:n] = (6.0 * (self.ydata[0:n - 1] - self.ydata[1:n]) /
                  (self.xdata[0:n - 1] - self.xdata[1:n])
                  - 6.0 * (self.ydata[1:n] - self.ydata[2:n + 1]) /
                  (self.xdata[1:n] - self.xdata[2:n + 1]))
        self.ludecomp3(c, d, e)
        self.lusolve3(c, d, e, k)
        return k

    def find_segment(self, x):
        iLeft = 0
        iRight = len(self.xdata) - 1
        while 1:
            if (iRight - iLeft) <= 1:
                return iLeft
            i = (iLeft + iRight) / 2
            if x < self.xdata[i]:
                iRight = i
            else:
                iLeft = i

    def eval_spline(self, x):
        """
        y = evalSpline(xData, yData, k, x)
        Evaluates cubic spline at x. The curvatures k can be computed with the
        function "curvatures".
        """
        i = self.find_segment(x)
        h = self.xdata[i] - self.xdata[i + 1]
        y = (((x - self.xdata[i + 1])**3.0/h - (x - self.xdata[i + 1])*h)*self.k[i] / 6.0
             - ((x - self.xdata[i])**3.0/h - (x - self.xdata[i])*h)*self.k[i + 1] / 6.0
             + (self.ydata[i]*(x - self.xdata[i + 1]) - self.ydata[i + 1]*(x - self.xdata[i])) / h)
        return y

    @property
    def xdata(self):
        return self._xdata

    @property
    def ydata(self):
        return self._ydata

    @property
    def k(self):
        return self._k


def linear_extrapolation(x, x1, x2, y1, y2):
    return y1 + (x - x1)*((y2 - y1)/(x2-x1))