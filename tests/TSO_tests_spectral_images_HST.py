# -*- coding: utf-8 -*-
"""
CASCADe

Test reading HST spectroscopic images

@author: Jeroen Bouwman
"""
import numpy as np
from matplotlib import pyplot as plt
from astropy.visualization import quantity_support
from astropy.io import ascii
from astropy import units as u

import cascade
from cascade.utilities import spectres

# create transit spectoscopy object
tso = cascade.TSO.TSOSuite()

# reset
tso.execute("reset")

# initialization with ini files
path = cascade.initialize.default_initialization_path
tso = cascade.TSO.TSOSuite("initialize", "cascade_test_HST_cpm.ini",
                           "cascade_test_HST_object.ini",
                           "cascade_test_HST_data_spectral_images.ini",
                           path=path)

assert tso.cascade_parameters == cascade.initialize.cascade_configuration
assert tso.cascade_parameters.isInitialized is True

# load data
tso.execute("load_data")

image0 = tso.observation.dataset.data[:, :, 0].copy()
image0.set_fill_value(np.nan)
plt.imshow(image0.filled().value,
           origin='lower',
           cmap='hot',
           interpolation='nearest',
           aspect='equal')
plt.colorbar().set_label("Intensity ({})".format(image0.data.unit))
plt.xlabel("Pixel number spatial direction")
plt.ylabel("Pixel number")
plt.title('Science data')
plt.show()

image1 = tso.observation.dataset_background.data[:, :, 0].copy()
image1.set_fill_value(np.nan)
plt.imshow(image1.filled().value,
           origin='lower',
           cmap='hot',
           interpolation='nearest',
           aspect='equal')
plt.colorbar().set_label("Intensity ({})".format(image1.data.unit))
plt.xlabel("Pixel number spatial direction")
plt.ylabel("Pixel number")
plt.title('Background data')
plt.show()

assert tso.observation.dataset.time.data.shape == (128, 128, 274)
assert type(tso.observation.dataset.time.data) == u.quantity.Quantity
assert tso.observation.dataset.time.data.unit == u.dimensionless_unscaled

# subtract background
tso.execute("subtract_background")

# sigma clip data
tso.execute("sigma_clip_data")

# eclipse model
tso.execute("define_eclipse_model")

# determine position of source from data set
tso.execute("determine_source_position")
with quantity_support():
    plt.plot(tso.observation.spectral_trace['wavelength_pixel'],
             tso.observation.spectral_trace['positional_pixel'] -
             np.median(tso.observation.spectral_trace['positional_pixel']))
    plt.plot(tso.cpm.spectral_trace)
    plt.show()

with quantity_support():
    plt.plot(tso.observation.spectral_trace['wavelength_pixel'],
             tso.observation.spectral_trace['wavelength'])
    plt.show()

plt.plot(tso.cpm.spectral_trace)
plt.show()

plt.plot(tso.observation.dataset.time.data.value[80, 85, :],
         tso.cpm.position[80, 85, :])
plt.show()

from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans
from skimage.feature import register_translation
image0 = np.ma.array(tso.observation.dataset.data.data.value[:, :, 0],
                     mask=tso.observation.dataset.mask[:, :, 0])
np.ma.set_fill_value(image0, float("NaN"))
kernel = Gaussian2DKernel(x_stddev=0.2, y_stddev=3.0, theta=-0.0)
cleaned_image0 = interpolate_replace_nans(image0.filled(), kernel)
plt.imshow(image0)
plt.show()
plt.imshow(cleaned_image0)
plt.show()

_, _, nintegrations = tso.observation.dataset.data.data.value.shape
shift_store = np.zeros((2, nintegrations))
for it in range(nintegrations):
    # subpixel precision
    image = np.ma.array(tso.observation.dataset.data.data.value[:, :, it],
                        mask=tso.observation.dataset.mask[:, :, it])
    np.ma.set_fill_value(image, float("NaN"))
    cleaned_image = interpolate_replace_nans(image.filled(), kernel)
    shift, error, diffphase = register_translation(cleaned_image0, cleaned_image, 150)
    shift_store[:, it] = shift

fig = plt.figure(figsize=(8, 4))
ax1 = plt.subplot(1, 2, 1, adjustable='box')
ax2 = plt.subplot(1, 2, 2, sharex=ax1, sharey=ax1, adjustable='box')
ax1.plot(tso.observation.dataset.time.data.value[80, 18, :], shift_store[0, :])
ax1.plot(tso.observation.dataset.time.data.value[80, 18, :],
         tso.observation.instrument_calibration.relative_source_shift['disp_shift'])
ax1.set_title('Offset in dispersion direction')
ax2.plot(tso.observation.dataset.time.data.value[80, 18, :], shift_store[1, :])
ax2.plot(tso.observation.dataset.time.data.value[80, 18, :],
         tso.observation.instrument_calibration.relative_source_shift['cross_disp_shift'])
ax2.set_title('Offset cross dispertion direction')
plt.show()

fig = plt.figure(figsize=(8, 4))
ax1 = plt.subplot(1, 1, 1)
ax1.plot(tso.observation.dataset.time.data.value[80, 18, :],
         tso.cpm.position[80, 18, :])
# ax1.plot(tso.observation.dataset.time.data.value[80, 18, :],
#         shift_store[0, :] - np.median(shift_store[0, :]))
ax1.plot(tso.observation.dataset.time.data.value[80, 18, :],
         -shift_store[1, :] - np.median(-shift_store[1, :]))
ax1.set_ylim([-0.1, 0.1])
ax1.set_title('Comparisson between C.O.L. and Phase Correlation')
ax1.set_xlabel("Phase")
ax1.set_ylabel("Relative Offset [pix]")
plt.show()

# set the extraction area
tso.execute("set_extraction_mask")

# setup regressors
tso.execute("select_regressors")

# create calibrated time series and derive planetary signal
tso.execute("calibrate_timeseries")

# extract planetary signal
tso.execute("extract_spectrum")

max_lambda = np.max(tso.observation.dataset.wavelength.data.value)
min_lambda = np.min(tso.observation.dataset.wavelength.data.value)
extent = [0, 127, max_lambda, min_lambda]
delta_lam = max_lambda - min_lambda
fig, ax = plt.subplots(figsize=(6, 6))
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
cmap = plt.cm.gist_heat
cmap.set_bad('black', 1.)
ax.imshow(np.ma.abs(tso.exoplanet_spectrum.weighted_image),
          origin='lower', aspect='auto',
          cmap=cmap, interpolation='none', extent=extent, vmin=0, vmax=1000)
ax.set_aspect(128.0 / delta_lam)
ax.set_xlabel('Pixel Number Spatial Direction')
ax.set_ylabel('Wavelength')

path_old = '/home/bouwman/'
spec_instr_model_mandell = ascii.read(path_old+'WASP19b_mandell.txt',
                                      data_start=0)
wave_mandell = (spec_instr_model_mandell['col1'] * u.micron)
flux_mandell = (spec_instr_model_mandell['col2'] * u.percent)
error_mandell = (spec_instr_model_mandell['col3'] * u.percent)

spec_instr_model_huitson = ascii.read(path_old+'WASP19b_huitson.txt',
                                      data_start=0)
wave_huitson = (spec_instr_model_huitson['col1']*u.micron)
flux_huitson = (spec_instr_model_huitson['col2'])
error_huitson = (spec_instr_model_huitson['col3'])
rel_error = error_huitson/flux_huitson
# the factor 1.007 to correct for limb darkning
flux_huitson = (flux_huitson*1.007)**2 * 100 * u.percent
error_huitson = np.sqrt(2)*rel_error*flux_huitson

# tso.exoplanet_spectrum.spectrum.wavelength_unit = u.micrometer
mask_use = tso.exoplanet_spectrum.spectrum.wavelength.mask

rebinned_wave = \
    tso.exoplanet_spectrum.spectrum.wavelength.data.value[~mask_use][3:-3:6]
rebinned_spec, rebinned_error = \
    spectres(rebinned_wave,
             tso.exoplanet_spectrum.spectrum.wavelength.data.value[~mask_use],
             tso.exoplanet_spectrum.spectrum.data.data.value[~mask_use],
             spec_errs=tso.exoplanet_spectrum.
             spectrum.uncertainty.data.value[~mask_use])

fig, ax = plt.subplots(figsize=(9, 7))
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
ax.plot(rebinned_wave, rebinned_spec, color="black", lw=5, alpha=0.9, zorder=5)
ax.errorbar(rebinned_wave, rebinned_spec, yerr=rebinned_error,
            fmt=".k", color="black", lw=5,
            alpha=0.9, ecolor="black",
            markeredgecolor='black', fillstyle='full', markersize=10,
            markerfacecolor='black', zorder=5)
ax.plot(wave_mandell, flux_mandell, color="r", lw=5, alpha=0.9, zorder=4)
ax.errorbar(wave_mandell.value, flux_mandell.value, yerr=error_mandell.value,
            fmt=".k", color="r", lw=5,
            alpha=0.9, ecolor="r",
            markeredgecolor='r', fillstyle='full', markersize=10,
            markerfacecolor='r', zorder=4)
ax.plot(wave_huitson, flux_huitson, color="y", lw=5, alpha=0.9, zorder=3)
ax.errorbar(wave_huitson.value, flux_huitson.value, yerr=error_huitson.value,
            fmt=".k", color="y", lw=5,
            alpha=0.9, ecolor="y",
            markeredgecolor='y', fillstyle='full', markersize=10,
            markerfacecolor='y', zorder=3)
ax.plot(tso.exoplanet_spectrum.spectrum.wavelength,
        tso.exoplanet_spectrum.spectrum.data, lw=3, alpha=0.5, color='gray',
        zorder=2)
ax.errorbar(tso.exoplanet_spectrum.spectrum.wavelength.data.value,
            tso.exoplanet_spectrum.spectrum.data.data.value,
            yerr=tso.exoplanet_spectrum.spectrum.uncertainty.data.value,
            fmt=".k", color='gray', lw=3, alpha=0.5, ecolor='gray',
            markerfacecolor='gray',
            markeredgecolor='gray', fillstyle='full', markersize=10,
            zorder=2)
axes = plt.gca()
axes.set_xlim([1.1, 1.7])
axes.set_ylim([1.9, 2.1])
ax.set_ylabel('Transit depth')
ax.set_xlabel('Wavelength')
plt.show()

# save planetary signal
tso.execute("save_results")

# plot planetary signal
tso.execute("plot_results")
