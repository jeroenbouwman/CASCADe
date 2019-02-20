import unittest
import cascade
from cascade.TSO import TSOSuite


class TestCascade(unittest.TestCase):
    def setUp(self):
        self.tso = TSOSuite()

    def tearDown(self):
        del self.tso

    def test_basic_cascade(self):
        self.assertIsInstance(self.tso, cascade.TSO.TSOSuite)


if __name__ == '__main__':
    #  unittest.main()
    suite = unittest.TestLoader().loadTestsFromTestCase(TestCascade)
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite)
