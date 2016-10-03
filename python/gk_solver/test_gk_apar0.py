import numpy.testing as npt
import numpy as np
import unittest
from gk_apar0 import dispersion, det_deriv

class Test_gk_apar0(unittest.TestCase):
    def test_det_deriv(self):
        # a set of test parameters
        ti_te = [10., 10, 10, 1, 1, 1, 0.1, 0.1, 0.1]
        kperp_rhoi = [0.5, 0.3, 0.1, 10, 1, 0.1, 0.9, 0.5, 0.1]
        bi = [0.1, 0.3, 0.5, 0.7, 0.9, 1.5, 3, 7, 10]
        mi_me = 1836
        wbar = [0.1, 0.3, 0.5, 0.9, 1.8, 3.6, 5, 7.2, 10]
        epsilon  = 1e-3

        for i in range(len(ti_te)):
            deriv_val = det_deriv(ti_te[i], mi_me, bi[i], kperp_rhoi[i], wbar[i])
            deriv_approx = (dispersion(ti_te[i], mi_me, bi[i], kperp_rhoi[i],
                                       wbar[i]+epsilon) - 
                            dispersion(ti_te[i], mi_me, bi[i], kperp_rhoi[i], 
                                       wbar[i]-epsilon))/(2 * epsilon)
            npt.assert_allclose(deriv_val, deriv_approx, rtol=1e-2)

if __name__ == '__main__':
    unittest.main()
