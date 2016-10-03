import numpy as np
from scipy.special import wofz

def real_imag(val):
    """
    Return the list [real(val), imag(val)]

    """
    return [np.real(val), np.imag(val)]

def list2complex(x):
    """
    Convert a list of two numbers to a complex number:
    x -> x[0] + x[1] * j
    """
    return x[0] + x[1] * 1j


def xi_i(bi, w_bar):
    """
    calculate dimensionless parameter for ion
    xi_i = p / (-i k_par v_ti )
    """
    return w_bar / np.sqrt(bi)

def xi_e(bi, w_bar, ti_te, mi_me):
    """
    calculate dimensionless parameter for ion
    xi_e = p / (-i k_par v_te )
    """
    return w_bar * np.sqrt(ti_te / mi_me / bi) 

def zp(x):
    """
    plasma dispersion function
    using faddeeva function from scipy.

    """
    # return -2. * dawsn(x) + 1j * np.sqrt(np.pi) * np.exp(- x**2)

    sqrt_pi = 1.7724538509055159
    return 1j * sqrt_pi * wofz(x)