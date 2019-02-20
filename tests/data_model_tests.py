# -*- coding: utf-8 -*-
import cascade
from cascade.data_model import SpectralData
from cascade.data_model import SpectralDataTimeSeries
import unittest
import numpy as np
import astropy.units as u


class TestDataModel(unittest.TestCase):
    def setUp(self):
        # initialization with wavelength and flux array
        # without units
        self.wave = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
        self.flux = np.array([8.0, 8.0, 8.0, 8.0, 8.0, 8.0])
        # with units
        self.wave_units = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])*u.micron
        self.flux_units = np.array([8.0, 8.0, 8.0, 8.0, 8.0, 8.0])*u.Jy
        self.time_units = np.array([240000.0, 2400001.0, 2400002.0])*u.day
        # mask
        self.mask = np.array([False, False, False, True, False, True])
        # auxilery data
        self.position_units = np.array([-0.1, 0.0, 0.1])*u.pix
        self.temperature_units = np.array([5.0, 6.0, 7.0])*u.K
        # different flux unit
        self.flux_unit_change = u.erg/u.s/u.cm**2/u.Hz
        self.wavelength_unit_change = u.Hz
        # flux changed to different unit
        self.changed_flux = self.flux_units.to(self.flux_unit_change)
        self.changed_wavelength = \
            self.wave_units.to(self.wavelength_unit_change,
                               equivalencies=u.spectral())

    def tearDown(self):
        del self.wave
        del self.flux
        del self.wave_units
        del self.flux_units
        del self.time_units
        del self.mask
        del self.position_units
        del self.temperature_units
        del self.flux_unit_change
        del self.wavelength_unit_change
        del self.changed_flux
        del self.changed_wavelength

    def test_basic_Data(self):
        # initialization without data
        sd = SpectralData()
        self.assertIsInstance(sd, cascade.data_model.SpectralData)

        # initialization without units
        sd = SpectralData(wavelength=self.wave,
                          data=self.flux)
        self.assertTrue((sd.data == self.flux).data.value.all())
        self.assertEqual(sd.wavelength_unit, None)

        # initialization with units
        sd = SpectralData(wavelength=self.wave_units,
                          data=self.flux_units)
        self.assertTrue((sd.data == self.flux_units).data.value.all())
        self.assertEqual(sd.wavelength_unit, self.wave_units.unit)
        # changing units
        sd.data_unit = self.flux_unit_change
        self.assertTrue((sd.data == self.changed_flux).data.value.all())

        # initialization with units specifying different unit
        sd = SpectralData(wavelength=self.wave_units,
                          data=self.flux_units,
                          wavelength_unit=self.wavelength_unit_change,
                          data_unit=self.flux_unit_change)
        self.assertTrue((sd.data == self.changed_flux).data.value.all())
        self.assertTrue((sd.wavelength == self.changed_wavelength).
                        data.value.all())

        # initialization with masked data
        masked_wave = np.ma.array(self.wave_units, mask=self.mask)
        masked_flux = np.ma.array(self.flux_units, mask=self.mask)
        sd = SpectralData(wavelength=masked_wave,
                          data=masked_flux)
        self.assertTrue((sd.mask == self.mask).all())
        self.assertTrue((sd.data.mask == self.mask).all())

    def test_basic_DataTimeSeries(self):
        # test spectral data time series class
        # besic init
        sdt = SpectralDataTimeSeries()
        self.assertIsInstance(sdt, cascade.data_model.SpectralDataTimeSeries)

        # initialization with units
        flux = np.repeat(self.flux_units[:, np.newaxis],
                         self.time_units.shape[0], 1)
        sdt = SpectralDataTimeSeries(wavelength=self.wave_units,
                                     data=flux,
                                     time=self.time_units)
        self.assertTrue((sdt.data == flux).data.value.all())
        self.assertEqual(sdt.time_unit, self.time_units.unit)

        # with auxilary data such as position or temperature
        sdt = SpectralDataTimeSeries(wavelength=self.wave_units,
                                     data=flux,
                                     time=self.time_units,
                                     position=self.position_units,
                                     temperature=self.temperature_units)
        self.assertTrue(sdt.position.ndim == 2)
        self.assertTrue((sdt.position == self.position_units).data.value.all())
        self.assertTrue(sdt.temperature.ndim == 2)
        self.assertTrue((sdt.temperature == self.temperature_units).
                        data.value.all())


if __name__ == '__main__':
    # unittest.main()
    suite = unittest.TestLoader().loadTestsFromTestCase(TestDataModel)
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite)
