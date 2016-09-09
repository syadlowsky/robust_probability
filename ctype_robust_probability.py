# coding: utf-8
import ctypes
from numpy.ctypeslib import ndpointer
import numpy as np

minimize_linear_square_divergence = ctypes.cdll.LoadLibrary(
    './c_robust_probability/robust_probability.so').minimize_linear_square_divergence

minimize_linear_square_divergence.argtypes = \
    [ctypes.c_double, ndpointer(ctypes.c_double),
     ctypes.c_int, ctypes.c_double, ctypes.c_double,
     ndpointer(ctypes.c_double)]
minimize_linear_square_divergence.restype = ctypes.c_double

def robust_probability(z, rho, scale=None, mass=1.0):
    scale = np.float64(scale or z.size)
    mass = np.float64(mass)
    out = np.zeros(z.size, dtype=np.float64)
    minimize_linear_square_divergence(rho, -z, z.size, mass, scale, out)
    return out
