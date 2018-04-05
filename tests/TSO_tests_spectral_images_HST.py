# -*- coding: utf-8 -*-

import cascade
from matplotlib import pyplot as plt
from astropy.visualization import quantity_support
import numpy as np
from astropy.io import ascii
from astropy import units as u
# from astropy.units import cds
from scipy import stats


# create transit spectoscopy object
tso = cascade.TSO.TSOSuite()

# reset
tso.execute("reset")

# initialization with ini files
path = cascade.initialize.default_initialization_path
tso = cascade.TSO.TSOSuite("initialize", "cascade_test_HST_cpm.ini",
                           "cascade_test_HST_object.ini",
                           "cascade_test_HST_data_spectral_images.ini", path=path)

print(tso.cascade_parameters)
print(cascade.initialize.cascade_configuration)
print(tso.cascade_parameters.isInitialized)

# load data
tso.execute("load_data")
plt.imshow(tso.observation.dataset.data.data.value[:, :, 0])
plt.show()
with quantity_support():
    plt.plot(tso.observation.dataset.time.data[80, 85, :],
             tso.observation.dataset.data.data[80, 85, :])
    plt.show()

plt.imshow(tso.observation.dataset.wavelength.data.value[:, :, 0])
plt.show()

# subtract background
tso.execute("subtract_background")
plt.imshow(tso.observation.dataset.data.data.value[:, :, 0])
plt.show()
with quantity_support():
    plt.plot(tso.observation.dataset.time[80, 85, :],
             tso.observation.dataset.data[80, 85, :])
    plt.show()

background_model_parameters = \
    tso.observation.instrument_calibration.background_model_parameters
plt.plot(background_model_parameters['parameter'])
plt.errorbar(range(background_model_parameters['parameter'].size),
             background_model_parameters['parameter'],
             yerr=background_model_parameters['error'])
plt.show()

# sigma clip data
tso.execute("sigma_clip_data")
plt.imshow(np.ma.median(tso.observation.dataset.data[:, :, :], axis=2))
plt.show()
plt.imshow(tso.observation.dataset.mask[:, :, :].all(axis=2))
plt.show()
with quantity_support():
    plt.plot(tso.observation.dataset.time[80, 85, :],
             tso.observation.dataset.data[80, 85, :])
    plt.show()

# eclipse model
tso.execute("define_eclipse_model")
plt.imshow(tso.model.light_curve_interpolated[0][:, 0, :])
plt.show()
with quantity_support():
    plt.plot(tso.observation.dataset.time.data[80, 85, :],
             tso.model.light_curve_interpolated[0][80, 85, :])
    plt.show()


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
                     mask = tso.observation.dataset.mask[:, :, 0])
np.ma.set_fill_value(image0, float("NaN"))
kernel = Gaussian2DKernel(x_stddev=1.5)
cleaned_image0 = interpolate_replace_nans(image0.filled(),
                                         kernel)
plt.imshow(image0)
plt.show()
plt.imshow(cleaned_image0)
plt.show()

_,_,nintegrations = tso.observation.dataset.data.data.value.shape
shift_store = np.zeros((2, nintegrations))
for it in range(nintegrations):
    # subpixel precision
    image = np.ma.array(tso.observation.dataset.data.data.value[:, :, it],
                     mask = tso.observation.dataset.mask[:, :, it])
    np.ma.set_fill_value(image, float("NaN"))
    cleaned_image = interpolate_replace_nans(image.filled(),
                                         kernel)
    shift, error, diffphase = register_translation(cleaned_image0, cleaned_image, 150)
    shift_store[:,it] = shift

fig = plt.figure(figsize=(8, 3))
ax1 = plt.subplot(1, 2, 1, adjustable='box-forced')
ax2 = plt.subplot(1, 2, 2, sharex=ax1, sharey=ax1, adjustable='box-forced')
ax1.plot(tso.observation.dataset.time.data.value[80, 18, :], shift_store[0,:])
ax1.set_title('Y offset')
ax2.plot(tso.observation.dataset.time.data.value[80, 18, :], shift_store[1,:])
ax2.set_title('X offset')
plt.show()

fig = plt.figure(figsize=(8, 3))
ax1 = plt.subplot(1, 1, 1)
ax1.plot(tso.observation.dataset.time.data.value[80, 18, :],
         tso.cpm.position[80, 18, :])
ax1.plot(tso.observation.dataset.time.data.value[80, 18, :],shift_store[0,:] - np.median(shift_store[0,:]))
ax1.plot(tso.observation.dataset.time.data.value[80, 18, :],-shift_store[1,:]- np.median(-shift_store[1,:]))
ax1.set_ylim([-0.1,0.1])
plt.show()

# set the extraction area
tso.execute("set_extraction_mask")

print(tso.cpm.extraction_mask[0].shape)
plt.imshow(tso.cpm.extraction_mask[0])
plt.show()

# setup regressors
tso.execute("select_regressors")


# create calibrated time series and derive planetary signal
tso.execute("calibrate_timeseries")
print(np.ma.max(tso.calibration_results.regularization))
print(np.ma.min(tso.calibration_results.regularization))
plt.imshow(tso.calibration_results.regularization)
plt.show()


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

path_old = '/home/bouwman//'
spec_instr_model_mandell = ascii.read(path_old+'WASP19b_mandell.txt', data_start=0)
wave_mandell = (spec_instr_model_mandell['col1']*u.micron)
flux_mandell = (spec_instr_model_mandell['col2'] *u.percent)
error_mandell = (spec_instr_model_mandell['col3'] * u.percent)

tso.exoplanet_spectrum.spectrum.wavelength_unit = u.micrometer

from pysynphot import observation
from pysynphot import spectrum

def rebin_spec(wave, specin, wavnew):
    spec = spectrum.ArraySourceSpectrum(wave=wave, flux=specin)
    f = np.ones(len(wave))
    filt = spectrum.ArraySpectralElement(wave, f, waveunits='micron')
    obs = observation.Observation(spec, filt, binset=wavnew, force='taper')
    print(obs)
    return obs.binflux

mask_use = tso.exoplanet_spectrum.spectrum.wavelength.mask
rebinned_spec = rebin_spec(tso.exoplanet_spectrum.spectrum.wavelength.data.value[~mask_use],
                           tso.exoplanet_spectrum.spectrum.data.data.value[~mask_use],
                           tso.exoplanet_spectrum.spectrum.wavelength.data.value[~mask_use][::6])
rebinned_wave = tso.exoplanet_spectrum.spectrum.wavelength.data.value[~mask_use][::6]
fig, ax = plt.subplots(figsize=(7, 4))
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
ax.plot(rebinned_wave, rebinned_spec, color="g", lw=5, alpha=0.9)
ax.plot(wave_mandell, flux_mandell, color="r", lw=3, alpha=0.9)
ax.errorbar(wave_mandell.value, flux_mandell.value, yerr=error_mandell.value,
            fmt=".k", color="r", lw=2,
            alpha=0.9, ecolor="r",
            markeredgecolor='r', fillstyle='full', markersize=10,
            markerfacecolor='r', zorder=3)
ax.plot(tso.exoplanet_spectrum.spectrum.wavelength,
        tso.exoplanet_spectrum.spectrum.data, lw=3, alpha=0.5, color='blue')
ax.errorbar(tso.exoplanet_spectrum.spectrum.wavelength.data.value,
            tso.exoplanet_spectrum.spectrum.data.data.value,
            yerr=tso.exoplanet_spectrum.spectrum.uncertainty.data.value,
            fmt=".k", color='blue', lw=3, alpha=0.5, ecolor='blue',
            markerfacecolor='blue',
            markeredgecolor='blue', fillstyle='full', markersize=10,
            zorder=4)
axes = plt.gca()
axes.set_xlim([1.1, 1.7])
axes.set_ylim([1.9, 2.2])
ax.set_ylabel('Transit depth')
ax.set_xlabel('Wavelength')
plt.show()


# save planetary signal
tso.execute("save_results")

# plot planetary signal
tso.execute("plot_results")

cal_signal_depth = float(tso.cascade_parameters.cpm_calibration_signal_depth)
mean_eclipse_depth = float(tso.cascade_parameters.observations_median_signal)

plt.plot(tso.calibration_results.calibrated_time_series[50, :].T)
plt.show()

plt.plot(tso.exoplanet_spectrum.spectrum.wavelength,
         tso.calibration_results.calibration_signal *
         (1+cal_signal_depth)+cal_signal_depth)
axes = plt.gca()
axes.set_ylim([cal_signal_depth*0.0, cal_signal_depth*5])
plt.show()

plt.plot(tso.exoplanet_spectrum.spectrum.wavelength,
         (tso.calibration_results.calibration_signal *
          (1+cal_signal_depth)+cal_signal_depth)/cal_signal_depth)
axes = plt.gca()
axes.set_ylim([0.0, 2.5])
plt.show()

plt.imshow(np.log10((tso.calibration_results.calibration_signal *
           (1+cal_signal_depth)+cal_signal_depth) /
           (tso.calibration_results.error_calibration_signal *
            (1+cal_signal_depth)/cal_signal_depth)**2))
plt.show()

plt.plot(tso.exoplanet_spectrum.wavelength,
         tso.exoplanet_spectrum.calibration_correction)
plt.errorbar(tso.exoplanet_spectrum.wavelength,
             tso.exoplanet_spectrum.calibration_correction,
             yerr=tso.exoplanet_spectrum.calibration_correction_error)
plt.show()

plt.plot(tso.exoplanet_spectrum.wavelength,
         (mean_eclipse_depth-tso.exoplanet_spectrum.calibration_correction) /
         tso.exoplanet_spectrum.calibration_correction_error)
plt.plot(tso.exoplanet_spectrum.wavelength,
         tso.exoplanet_spectrum.data/tso.exoplanet_spectrum.error)
plt.show()
