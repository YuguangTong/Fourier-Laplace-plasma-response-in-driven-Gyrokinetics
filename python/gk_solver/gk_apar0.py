# package to solve gyrokinetic equations when Apar = 0

import numpy as np
from scipy.special import i0e, i1e
import scipy.optimize
from gk_solver.util import real_imag, list2complex, zp

#-------------------------------------
# Terms ABCDE in the dispersion tensor
#-------------------------------------
def A(ti_te, mi_me, bi, kperp_rhoi, w_bar):
    """
    Calculate
    A  = sum_s [ (Ti/Ts) (1 + Gamma_0s xi_s Z_s) ]
    """
    alpha_i = kperp_rhoi**2 / 2
    alpha_e = alpha_i / ti_te / mi_me
    xi_i = w_bar / np.sqrt(bi)
    xi_e = w_bar * np.sqrt(ti_te / mi_me / bi)
    Z_i = zp(xi_i)
    Z_e = zp(xi_e)
    Gamma_0i = i0e(alpha_i)
    Gamma_0e = i0e(alpha_e)
    return 1 + Gamma_0i * xi_i * Z_i + ti_te * (1 + Gamma_0e * xi_e * Z_e)

def B(ti_te, mi_me, kperp_rhoi):
    """
    Calculate
    B = sum_s [ (Ti/Ts) ( 1 - Gamma_0(alpha_s)) ]
    """
    alpha_i = kperp_rhoi**2 / 2
    alpha_e = alpha_i / ti_te / mi_me
    return 1 - i0e(alpha_i) - ti_te * (1 - i0e(alpha_e))

def C(ti_te, mi_me, bi, kperp_rhoi, w_bar):
    """
    C = sum_s (qi/qs) Gamma_1s xi_s Z_s 

    """
    alpha_i = kperp_rhoi**2 / 2
    alpha_e = alpha_i / ti_te / mi_me
    xi_i = w_bar / np.sqrt(bi)
    xi_e = w_bar * np.sqrt(ti_te / mi_me / bi)
    Z_i = zp(xi_i)
    Z_e = zp(xi_e)
    Gamma_1i = i0e(alpha_i) - i1e(alpha_i)
    Gamma_1e = i0e(alpha_e) - i1e(alpha_e)
    
    res = Gamma_1i * xi_i * Z_i - Gamma_1e * xi_e * Z_e
    return res

def D(ti_te, mi_me, bi, kperp_rhoi, w_bar):
    """
    D = sum_s (2Ts/Ti) Gamma_1s xi_s Z_s 

    """
    alpha_i = kperp_rhoi**2 / 2
    alpha_e = alpha_i / ti_te / mi_me
    xi_i = w_bar / np.sqrt(bi)
    xi_e = w_bar * np.sqrt(ti_te / mi_me / bi)
    Z_i = zp(xi_i)
    Z_e = zp(xi_e)
    Gamma_1i = i0e(alpha_i) - i1e(alpha_i)
    Gamma_1e = i0e(alpha_e) - i1e(alpha_e)
    
    res = Gamma_1i * xi_i * Z_i + 1/ti_te * Gamma_1e * xi_e * Z_e
    res *= 2
    return res

def E(ti_te, mi_me, kperp_rhoi):
    """
    E = sum_s (qs/qi) * Gamma_1s
    """
    alpha_i = kperp_rhoi**2 / 2
    alpha_e = alpha_i / ti_te / mi_me    
    Gamma_1i = i0e(alpha_i) - i1e(alpha_i)
    Gamma_1e = i0e(alpha_e) - i1e(alpha_e)   
    return Gamma_1i - Gamma_1e

def F(ti_te, mi_me, kperp_rhoi):
    """
    F = sum_s (2Ts/Ti) * Gamma_1s
    """
    alpha_i = kperp_rhoi**2 / 2
    alpha_e = alpha_i / ti_te / mi_me    
    Gamma_1i = i0e(alpha_i) - i1e(alpha_i)
    Gamma_1e = i0e(alpha_e) - i1e(alpha_e)   
    return 2 * Gamma_1i - 2 * Gamma_1e / ti_te

def zp_prime(xi):
    """
    zp_prime = -2 * (1 + xi * zp(xi))
    """
    return -2 * (1 + xi * zp(xi))
    
def L(ti_te, mi_me, bi, w_bar):
    """
    L = sum_s i (qs/qi) *  (Ti/Ts) * (ms/mi) * Zp_prime(xs)
    """
    xi_i = w_bar / np.sqrt(bi)
    xi_e = w_bar * np.sqrt(ti_te / mi_me / bi)
    Zprime_i = zp_prime(xi_i)
    Zprime_e = zp_prime(xi_e)
    return 1j * (Zprime_i - ti_te/mi_me * Zprime_e)

def M(ti_te, mi_me, bi, w_bar):
    """
    M = sum_s i ms/mi Zp_prime(xs)
    """
    xi_i = w_bar / np.sqrt(bi)
    xi_e = w_bar * np.sqrt(ti_te / mi_me / bi)
    Zprime_i = zp_prime(xi_i)
    Zprime_e = zp_prime(xi_e)
    
    return 1j * (Zprime_i + 1/mi_me * Zprime_e)

#--------------------
# Dispersion relation
#--------------------
# Use Poisson equation and the perpendicular component
# of the Ampere's law

def dispersion(ti_te, mi_me, bi, kperp_rhoi, w_bar):
    """
    dispersion = det(disp tensor) = A*(D-2/beta_i) - C^2
    """
    a = A(ti_te, mi_me, bi, kperp_rhoi, w_bar)
    c = C(ti_te, mi_me, bi, kperp_rhoi, w_bar)
    d = D(ti_te, mi_me, bi, kperp_rhoi, w_bar)
    return a * (d -2/bi) - c**2


#----------------------------------
# Derivative of dispersion relation
#----------------------------------

# d(DetM)/d(\overline{\omega})
def det_deriv(ti_te, mi_me, bi, kperp_rhoi, w_bar):
    """
    evaluate derivative of Det M w.r.t. omega_bar
    where omega_bar = w/k_par/v_A
    """
    alpha_i = kperp_rhoi**2 / 2
    alpha_e = alpha_i / ti_te / mi_me
    xi_i = w_bar / np.sqrt(bi)
    xi_e = w_bar * np.sqrt(ti_te / mi_me / bi)
    Z_i = zp(xi_i)
    Z_e = zp(xi_e)
    Gamma_1i = i0e(alpha_i) - i1e(alpha_i)
    Gamma_1e = i0e(alpha_e) - i1e(alpha_e)
    Gamma_0i = i0e(alpha_i)
    Gamma_0e = i0e(alpha_e)
    
    G_i = 1/np.sqrt(bi) * ((1-2*xi_i**2)*Z_i - 2 * xi_i)
    G_e = 1/np.sqrt(bi) * ((1-2*xi_e**2)*Z_e - 2 * xi_e) \
        * np.sqrt(ti_te / mi_me)
    
    # Aprime, dA / d omega_bar
    Ap = Gamma_0i * G_i + ti_te * Gamma_0e * G_e
    Cp = Gamma_1i * G_i - Gamma_1e * G_e
    Dp = 2 * (Gamma_1i  * G_i + 1/ti_te * Gamma_1e * G_e)

    a = A(ti_te, mi_me, bi, kperp_rhoi, w_bar)
    c = C(ti_te, mi_me, bi, kperp_rhoi, w_bar)
    d = D(ti_te, mi_me, bi, kperp_rhoi, w_bar)    
    
    res = Ap * (d - 2/bi) + a * Dp - 2 * c * Cp
    return res


#---------
# Residues
#---------
# Residues of the Bromwich integral
# Arise from inverse Laplace transform of the 
# Fourier-Laplace solution

def res_i(ti_te, mi_me, bi, kperp_rhoi, wbar_0, wbar_i, tbar, turnoff=None):
    """
    Calculate Residue(p_i) or Residue(\overline{\omega}_i)
    where i >= 1

    turnoff: controls which equations contain antenna driving terms
    """
    a = A(ti_te, mi_me, bi, kperp_rhoi, wbar_i)
    c = C(ti_te, mi_me, bi, kperp_rhoi, wbar_i)
    d = D(ti_te, mi_me, bi, kperp_rhoi, wbar_i)
    f = F(ti_te, mi_me, kperp_rhoi)
    e = E(ti_te, mi_me, kperp_rhoi)
    dmdw = det_deriv(ti_te, mi_me, bi, kperp_rhoi, wbar_i)
    
    if turnoff == None:
        numer = np.array([(d - 2/bi)*e - c*(f+2/bi),\
                              -c*e + a*(f+2/bi)])
    if turnoff == 'poisson':
        numer = np.array([ - c*(f+2/bi), a*(f+2/bi)])
    elif turnoff == 'ampere':
        numer = np.array([(d - 2/bi)*e, -c*e])
    denom = dmdw * (wbar_i - wbar_0)
    return numer / denom * np.exp(-1j * wbar_i * tbar)

def res_0(ti_te, mi_me, bi, kperp_rhoi, wbar_0, tbar, turnoff=None):
    """
    Calculate Residue(p_0) or Residue(\overline{\omega}_0)
    """
    a = A(ti_te, mi_me, bi, kperp_rhoi, wbar_0)
    c = C(ti_te, mi_me, bi, kperp_rhoi, wbar_0)
    d = D(ti_te, mi_me, bi, kperp_rhoi, wbar_0)
    f = F(ti_te, mi_me, kperp_rhoi)
    e = E(ti_te, mi_me, kperp_rhoi)
    detm = dispersion(ti_te, mi_me, bi, kperp_rhoi, wbar_0)

    if turnoff == None:
        numer = np.array([(d - 2/bi)*e - c*(f+2/bi),\
                              -c*e + a*(f+2/bi)])
    if turnoff == 'poisson':
        numer = np.array([ - c*(f+2/bi), a*(f+2/bi)])
    elif turnoff == 'ampere':
        numer = np.array([(d - 2/bi)*e, -c*e])

    return numer / detm * np.exp(-1j * wbar_0 * tbar)

def unittest():
    """
    run unittests
    """
    filename = 'test_gk_apar0.py'
    np.testing.run_module_suite(file_to_run=filename)

