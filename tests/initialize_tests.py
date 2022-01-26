import unittest
import os
from cascade.initialize import cascade_configuration
from cascade.initialize import generate_default_initialization
from cascade.initialize import cascade_default_initialization_path
from cascade.initialize import configurator
from cascade.TSO import TSOSuite


class TestInitialize(unittest.TestCase):
    def setUp(self):
        self.tso = TSOSuite()
        self.path = cascade_default_initialization_path

    def tearDown(self):
        del self.tso
        del self.path

    def test_basic_cascade(self):
        # generate default initialization file and check if filr is generated
        generate_default_initialization()
        self.assertTrue(os.path.exists(self.path / "cascade_default.ini"))
        # initialize cascade unsing the default ini file and check if
        # initilization is successful. Also check if the instance of the
        # configurator is identical to that of the calss instance (singleton)
        cascade_param_configuration = \
            configurator(self.path / "cascade_default.ini")
        self.assertTrue(cascade_param_configuration.isInitialized)
        self.assertEqual(cascade_param_configuration, cascade_configuration)
        # check if instances of configurator are equal also in tso object
        self.assertEqual(self.tso.cascade_parameters, cascade_configuration)

        # reset initialiszation
        cascade_param_configuration.reset()
        self.assertFalse(cascade_param_configuration.isInitialized)
        self.assertFalse(self.tso.cascade_parameters.isInitialized)

        # re-initialize tso object and check if instaces are stil equal
        self.tso.execute("initialize", "cascade_default.ini", path=self.path)
        self.assertTrue(self.tso.cascade_parameters.isInitialized)
        self.assertEqual(self.tso.cascade_parameters, cascade_configuration)


if __name__ == '__main__':
    #  unittest.main()
    suite = unittest.TestLoader().loadTestsFromTestCase(TestInitialize)
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite)
