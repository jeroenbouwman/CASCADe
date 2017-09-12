# -*- coding: utf-8 -*-
import cascade
import matplotlib.pyplot as plt
import numpy as np
from astropy.visualization import quantity_support
import astropy.units as u

wave = np.arange(4, 15, 0.05) * u.micron
temp = 300 * u.K
flux = cascade.exoplanet_tools.Planck(wave, temp)

with quantity_support():
    plt.plot(wave, flux)
    plt.show()

# test downloading catalog
ct = cascade.exoplanet_tools.parse_database('EXOPLANETS.ORG', update=True)
# test extracting data record for single system
dr = cascade.exoplanet_tools.extract_exoplanet_data(ct, 'HD 189733 b')

print(dr)

# test generation of ligthcurve model
# first generate standard ini file and initialize cascade
cascade.initialize.generate_default_initialization()
path = cascade.initialize.default_initialization_path
cascade_param = cascade.initialize.configurator(path+"cascade_default.ini")

# define ligthcurve model
lc_model = cascade.exoplanet_tools.lightcuve()
print(lc_model.valid_models)
print(lc_model.par)

fig, axs = plt.subplots(1, 1, figsize=(12, 10))
axs.plot(lc_model.lc[0], lc_model.lc[1])
axs.set_ylabel(r'Normalized Signal')
axs.set_xlabel(r'Phase')
axes = plt.gca()
axes.set_xlim([0, 1])
axes.set_ylim([-1.1, 0.1])
plt.show()

# calculate brightness temperature
wave = np.array([5.0, 6.0, 7.0, 8.0, 9.0])*u.micron
contrast = np.array([0.005, 0.005, 0.005, 0.005, 0.005])*u.dimensionless_unscaled
bt = cascade.exoplanet_tools.\
 convert_spectrum_to_brighness_temperature(wave,
                                           contrast,
                                           5000.0*u.K, 0.7*u.R_sun,
                                           1.2*u.R_jup)
print(bt)

wave = np.array([5.0, 6.0, 7.0, 8.0, 9.0])*u.micron
contrast = np.array([0.005, 0.005, 0.005, 0.005, 0.005])*u.dimensionless_unscaled
mask = np.array([False, False, True, False, True])
masked_wave = np.ma.array(wave, mask=mask)
masked_contrast = np.ma.array(contrast, mask=mask)
bt = cascade.exoplanet_tools.\
 convert_spectrum_to_brighness_temperature(masked_wave,
                                           masked_contrast,
                                           5000.0*u.K, 0.7*u.R_sun,
                                           1.2*u.R_jup)
print(bt)

wave = np.array([5.0, 6.0, 7.0, 8.0, 9.0])*u.micron
contrast = np.array([0.005, 0.005, 0.005, 0.005,
                     0.005])*u.dimensionless_unscaled
error_contrast = np.array([0.001, 0.001, 0.001, 0.001,
                           0.001])*u.dimensionless_unscaled
mask = np.array([False, False, True, False, True])
masked_wave = np.ma.array(wave, mask=mask)
masked_contrast = np.ma.array(contrast, mask=mask)
masked_error = np.ma.array(error_contrast, mask=mask)
bt = cascade.exoplanet_tools.\
 convert_spectrum_to_brighness_temperature(masked_wave,
                                           masked_contrast,
                                           5000.0*u.K, 0.7*u.R_sun,
                                           1.2*u.R_jup,
                                           error=masked_error)
print(bt)



