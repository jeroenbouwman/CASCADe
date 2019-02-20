# -*- coding: utf-8 -*-
import cascade
from cascade.cpm_model import solve_linear_equation
import unittest
import numpy as np
from scipy.stats import norm as norm_stats


class TestCpmModel(unittest.TestCase):
    def setUp(self):
        # create linear system with RV
        n = 1000
        x1 = norm_stats.rvs(0, 1, size=n)
        x2 = norm_stats.rvs(0, 1, size=n)
        x3 = norm_stats.rvs(0, 1, size=n)
        self.answer = np.array([10.0, 40.0, 0.1])
        self.A = np.column_stack([x1, x2, x3])
        self.b = self.answer[0] * x1 + self.answer[1] * x2 + \
            self.answer[2] * x3

    def tearDown(self):
        del self.answer
        del self.A
        del self.b

    def test_basic_cpm(self):
        regularization_pararameters = {"lam0": 1.e-16, "lam1": 1.e-8,
                                       "nlam": 150}
        # solve linear Eq.
        P, Perr, opt_reg_par, _, _, _ = \
            solve_linear_equation(self.A, self.b, cv_method='gcv',
                                  reg_par=regularization_pararameters)
        for i, (result, error) in enumerate(zip(P, Perr)):
            self.assertAlmostEqual(result, self.answer[i], places=None,
                                   msg=None, delta=100*error)
        self.assertLess(opt_reg_par, regularization_pararameters["lam1"])
        self.assertGreater(opt_reg_par, regularization_pararameters["lam0"])


if __name__ == '__main__':
    #  unittest.main()
    suite = unittest.TestLoader().loadTestsFromTestCase(TestCpmModel)
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite)
