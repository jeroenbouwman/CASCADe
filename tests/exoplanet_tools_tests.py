import unittest
import os
import numpy as np
import astropy.units as u
from astropy.table.table import Table
import collections
from urllib.error import URLError
from socket import gaierror
import warnings
import cascade
from cascade.exoplanet_tools import parse_database
from cascade.exoplanet_tools import extract_exoplanet_data
from cascade.exoplanet_tools import convert_spectrum_to_brighness_temperature
from cascade.exoplanet_tools import lightcurve
from cascade.initialize import cascade_default_initialization_path
from cascade.initialize import cascade_default_path
from cascade.initialize import generate_default_initialization
from cascade.initialize import configurator


class TestExoplanetsTools(unittest.TestCase):
    def setUp(self):
        self.wave = np.array([4, 15]) * u.micron
        self.BBTemp = 300 * u.K
        self.test_flux = np.array([3.85319124e-11, 5.01642712e-09])
        self.test_flux_unit = u.erg / (u.cm**2 * u.Hz * u.s * u.sr)
        self.contrast = np.array([0.005, 0.005])*u.dimensionless_unscaled
        self.error_contrast = np.array([0.001, 0.001])*u.dimensionless_unscaled
        self.mask = np.array([False, False])
        self.test_bt = np.array([1780.54871058, 1143.99801537])
        self.test_bt_unit = u.K
        self.path = cascade_default_path
        self.catalog_name = ['EXOPLANETS.ORG', "TEPCAT",
                             "NASAEXOPLANETARCHIVE", "EXOPLANETS_A"]
        self.catalog_dir = ["exoplanets.org", "tepcat",
                            "NASAEXOPLANETARCHIVE", "EXOPLANETS_A"]
        self.catalog_file_name = ["exoplanets.csv", "allplanets.csv",
                                  "nasaexoplanetarchive.csv",
                                  "exoplanets_a.csv"]
        self.test_system_name = ['HD 189733 b', 'HD_189733', 'HD 189733 b',
                                 'HD 189733 b']
        self.test_search_radius = 1.0*u.arcminute
        self.path_init_files = cascade_default_initialization_path

    def tearDown(self):
        del self.wave
        del self.BBTemp
        del self.test_flux
        del self.test_flux_unit
        del self.contrast
        del self.error_contrast
        del self.mask
        del self.path
        del self.catalog_name
        del self.catalog_dir
        del self.catalog_file_name
        del self.test_system_name
        del self.path_init_files

    def test_basic_cascade(self):
        # test BB function
        BBflux = cascade.exoplanet_tools.planck(self.wave, self.BBTemp)
        self.assertTrue(BBflux.unit == self.test_flux_unit)

        # test brightness temperature
        masked_wave = np.ma.array(self.wave, mask=self.mask)
        masked_contrast = np.ma.array(self.contrast, mask=self.mask)
        masked_error = np.ma.array(self.error_contrast, mask=self.mask)
        brightness_temperature = \
            convert_spectrum_to_brighness_temperature(masked_wave,
                                                      masked_contrast,
                                                      5000.0*u.K, 0.7*u.R_sun,
                                                      1.2*u.R_jup,
                                                      error=masked_error)
        self.assertIsInstance(brightness_temperature, tuple)
        self.assertIsInstance(brightness_temperature[0],
                              np.ma.core.MaskedArray)
        self.assertTrue(brightness_temperature[0].data.unit ==
                        self.test_bt_unit)
        for i in range(len(self.test_bt)):
            self.assertAlmostEqual(brightness_temperature[0].data.value[i],
                                   self.test_bt[i], places=3)

        # exoplanet catalog tests
        for (system_name, name,
             directory, file_name) in zip(self.test_system_name,
                                          self.catalog_name,
                                          self.catalog_dir,
                                          self.catalog_file_name):
            # test downloading catalog
            try:
                catalog = parse_database(name, update=True)
            except (URLError, gaierror):
                warnings.warn("{} not reachable, "
                              "skipping archive in test".format(name))
                continue
            self.assertTrue(os.path.exists(os.path.join(self.path,
                                                        "exoplanet_data",
                                                        directory,
                                                        file_name)))
            self.assertIsInstance(catalog, list)
            self.assertIsInstance(catalog[0], Table)
            self.assertTrue(len(catalog[0]) > 1)
            # test extracting data record for single system
            data_record = \
                extract_exoplanet_data(catalog, system_name,
                                       search_radius=self.test_search_radius)
            self.assertIsInstance(catalog, list)
            self.assertIsInstance(catalog[0], Table)
            self.assertTrue(data_record[0]['NAME'][0] == system_name)

        # test lightcurve model using parameters from default ini file
        generate_default_initialization()
        cascade_param = \
            configurator(self.path_init_files / "cascade_default.ini")
        cascade_param.telescope_collecting_area = '4.525 m2'
        cascade_param.instrument_dispersion_scale = '1000 Angstrom'
        lc_model = lightcurve(cascade_param)
        self.assertIsInstance(lc_model.par, collections.OrderedDict)
        self.assertTrue('batman' in lc_model._lightcurve__valid_models)
        self.assertIsInstance(lc_model.lc, tuple)


if __name__ == '__main__':
    #  unittest.main()
    suite = unittest.TestLoader().loadTestsFromTestCase(TestExoplanetsTools)
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite)
