# -*- coding: utf-8 -*-
import cascade
import numpy as np
import astropy.units as u
# from matplotlib import pyplot as plt

# test spectral data class
# basic creation
sd = cascade.data_model.SpectralData()

# initialization with wavelength and flux array
wave = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
flux = np.array([8.0, 8.0, 8.0, 8.0, 8.0, 8.0])
sd = cascade.data_model.SpectralData(wavelength=wave, data=flux)

# initialization with units
wave = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])*u.micron
flux = np.array([8.0, 8.0, 8.0, 8.0, 8.0, 8.0])*u.Jy
sd = cascade.data_model.SpectralData(wavelength=wave, data=flux)
print(sd.data)
print(sd.wavelength)
# changing units
sd.data_unit = u.erg/u.s/u.cm**2/u.Hz
print(sd.data)

# initialization with units specifying different unit
wave = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])*u.micron
flux = np.array([8.0, 8.0, 8.0, 8.0, 8.0, 8.0])*u.Jy
new_data_unit = u.erg/u.s/u.cm**2/u.Hz
new_wav_unit = u.Hz
sd = cascade.data_model.SpectralData(wavelength=wave, data=flux,
                                     wavelength_unit=new_wav_unit,
                                     data_unit=new_data_unit)
print(sd.data)

mask = np.array([False, False, False, True, False, True])
masked_wave = np.ma.array(wave, mask=mask)
masked_flux = np.ma.array(flux, mask=mask)
sd = cascade.data_model.SpectralData(wavelength=masked_wave,
                                     data=masked_flux)
print(sd.data)
mask[0] = True
sd = cascade.data_model.SpectralData(wavelength=masked_wave,
                                     data=masked_flux,
                                     mask = mask)
print(sd.data)

error = 0.1*flux
masked_error = np.ma.array(error, mask=mask)
sd = cascade.data_model.SpectralData(wavelength=masked_wave,
                                     data=masked_flux,
                                     uncertainty = error)
print(sd.data)
print(sd.data_unit)
print(sd.wavelength_unit)
print(sd.wavelength)
print(sd.uncertainty)


# with image data
intensity = np.zeros((10, 10))
intensity[:, 5] = 1.0
# plt.imshow(intensity)
wave_range = np.arange(10)
wave = np.repeat(wave_range[:, np.newaxis], wave_range.shape[0], 1)
sd = cascade.data_model.SpectralData(wavelength=wave, data=intensity)

# test spectral data time series class
# besic init
sdt = cascade.data_model.SpectralDataTimeSeries()

# initialization with units
wave = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])*u.micron
flux = np.array([8.0, 8.0, 8.0, 8.0, 8.0, 8.0])*u.Jy
time = np.array([240000.0, 2400001.0, 2400002.0])*u.day
flux = np.repeat(flux[:, np.newaxis], time.shape[0], 1)
sdt = cascade.data_model.SpectralDataTimeSeries(wavelength=wave, data=flux,
                                                time=time)

# initialization with spectral images
time = np.array([240000.0, 2400001.0, 2400002.0])*u.day
intensity = np.zeros((10, 10))
intensity[:, 5] = 1.0
intensity_tso = np.dstack([intensity]*time.shape[0])*u.adu/u.s/u.pix
wave_range = np.arange(10)
wave = np.repeat(wave_range[:, np.newaxis], wave_range.shape[0], 1)*u.micron
simt = cascade.data_model.SpectralDataTimeSeries(wavelength=wave,
                                                 data=intensity_tso, time=time)
print(simt.data)

# with auxilary data such as position or temperature
position = np.ones(time.shape[0])
position_unit = u.pix
temperature = np.ones(time.shape[0])*6.0 * u.K
simt = cascade.data_model.SpectralDataTimeSeries(wavelength=wave,
                                                 data=intensity_tso, time=time,
                                                 position=position,
                                                 position_unit=position_unit,
                                                 temperature=temperature)
print(simt.time)
print(simt.temperature)
print(simt.position)
