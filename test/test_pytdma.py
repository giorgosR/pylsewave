import unittest
import numpy as np
from pylsewave.cynum import pytdma, cpwgrad


class TestPytdmaFunction(unittest.TestCase):
    def setUp(self):
        self.func = pytdma

    def test_1_tdma(self):
        A = np.array([[10, 2, 0, 0], [3, 10, 4, 0], [0, 1, 7, 5], [0, 0, 3, 4]],
                     dtype=float)
        a = np.array([3., 1, 3])
        b = np.array([10., 10., 7., 4.])
        c = np.array([2., 4., 5.])
        d = np.array([3, 4, 5, 6.])
        x = np.zeros(d.shape[0])
        self.func(a, b, c, d, x)
        return np.testing.assert_allclose(x, np.linalg.solve(A, d), rtol=1e-9, atol=0)


    def test_2_tdma(self):
        A = np.array([[10, 2, 0, 0], [3, 10, 4, 0], [0, 1, 7, 5], [0, 0, 3, 4]],
                     dtype=float)
        a = np.array([3., 1, 3])
        b = np.array([10., 10., 7., 4.])
        c = np.array([2., 4., 5.])
        d = np.array([3, 4, 5, 6.])
        x = np.zeros(d.shape[0])
        self.assertEqual(self.func(a, b, c, d, x), 0)


class TestPygradFunction(unittest.TestCase):
    def setUp(self):
        self.func = cpwgrad

    def test_1_grad(self):
        f = np.array([1, 2, 4, 7, 11, 16], dtype=float)
        x = np.zeros(f.shape[0])
        dx = 2.0
        f_x = np.gradient(f, dx)
        self.func(f, x, dx)
        return np.testing.assert_allclose(x, f_x, rtol=1e-12, atol=0)

    def test_2_grad(self):
        f = np.array([1, 2, 4, 7, 11, 16], dtype=float)
        x = np.zeros(f.shape[0])
        dx = 2.0
        self.assertEqual(self.func(f, x, dx), 0)

    def test_3_grad(self):
        f = np.array([1, 2, 4, 7, 11, 16], dtype=float)
        x = np.zeros(f.shape[0])
        dx = 2.0
        self.assertIsInstance(self.func(f, x, dx), int)


if __name__ == "__main__":
    unittest.main()
