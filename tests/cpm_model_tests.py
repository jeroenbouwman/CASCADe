# -*- coding: utf-8 -*-
import cascade
from cascade.cpm_model import ols
from cascade.cpm_model import return_lambda_grid
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
        # solve linear Eq.
        (P, Perr, _) = ols(self.A, self.b)
        for i, (result, error) in enumerate(zip(P, Perr)):
            self.assertAlmostEqual(result, self.answer[i], places=None,
                                   msg=None, delta=1.e4*error)


if __name__ == '__main__':
    #  unittest.main()
    suite = unittest.TestLoader().loadTestsFromTestCase(TestCpmModel)
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite)
