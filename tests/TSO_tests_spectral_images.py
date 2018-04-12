# -*- coding: utf-8 -*-
"""
CASCADe

Test reading Spitzer spectroscopic images

@author: Jeroen Bouwman
"""
import cascade
from matplotlib import pyplot as plt
from astropy.visualization import quantity_support
import numpy as np
from astropy.io import ascii
from astropy import units as u
from scipy import stats

# create transit spectoscopy object
tso = cascade.TSO.TSOSuite()

# initialization with ini files
path = cascade.initialize.default_initialization_path
tso = cascade.TSO.TSOSuite("cascade_test_cpm.ini",
                           "cascade_test_object.ini",
                           "cascade_test_data_spectral_images.ini", path=path)
print(tso.cascade_parameters)
print(cascade.initialize.cascade_configuration)
print(tso.cascade_parameters.isInitialized)
assert tso.cascade_parameters == cascade.initialize.cascade_configuration
assert tso.cascade_parameters.isInitialized is True

# reset parameters
tso.execute("reset")
print(tso.cascade_parameters.isInitialized)
assert tso.cascade_parameters.isInitialized is False

# initialize without providing ini files
tso.execute("initialize")
print(tso.cascade_parameters.isInitialized)
assert tso.cascade_parameters.isInitialized is False

# create TSO object and initialize providing ini files
tso = cascade.TSO.TSOSuite()
path = cascade.initialize.default_initialization_path
tso.execute("initialize", "cascade_test_cpm.ini", "cascade_test_object.ini",
            "cascade_test_data_spectral_images.ini", path=path)
print(tso.cascade_parameters)
print(cascade.initialize.cascade_configuration)
print(tso.cascade_parameters.isInitialized)
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

image0 = tso.observation.dataset._wavelength.copy()
plt.imshow(image0,
           origin='lower',
           cmap='hot',
           interpolation='nearest',
           aspect='equal')
plt.colorbar().set_label("Wavelength ({})".
                         format(tso.observation.dataset.wavelength_unit))
plt.xlabel("Pixel number spatial direction")
plt.ylabel("Pixel number")
plt.title('Wavelength')
plt.show()

image1 = tso.observation.dataset.wavelength[:, :, 0]
image1.set_fill_value(np.nan)
plt.imshow(image1.filled().value,
           origin='lower',
           cmap='hot',
           interpolation='nearest',
           aspect='equal')
plt.colorbar().set_label("Intensity ({})".format(image1.data.unit))
plt.xlabel("Pixel number spatial direction")
plt.ylabel("Pixel number")
plt.title('Wavelength')
plt.show()

with quantity_support():
    plt.plot(tso.observation.dataset.time.data[80, 18, :],
             tso.observation.dataset.data.data[80, 18, :])
    plt.show()

with quantity_support():
    plt.plot(tso.observation.dataset.wavelength.data[:, 18, 0],
             tso.observation.dataset.data.data[:, 18, 0])
    plt.show()

print(tso.observation.dataset.time.data.shape)
print(type(tso.observation.dataset.time.data))
print(tso.observation.dataset.time.data.unit)
assert tso.observation.dataset.time.data.shape == (128, 128, 280)
assert type(tso.observation.dataset.time.data) == u.quantity.Quantity
assert tso.observation.dataset.time.data.unit == u.dimensionless_unscaled

# subtract background
tso.execute("subtract_background")

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

plt.imshow(tso.observation.dataset.mask.any(axis=2),
           origin='lower',
           cmap='hot',
           interpolation='nearest',
           aspect='equal')
plt.xlabel("Pixel number spatial direction")
plt.ylabel("Pixel number")
plt.title('Collapsed mask science data')
plt.show()

with quantity_support():
    plt.plot(tso.observation.dataset.time.data[80, 18, :],
             tso.observation.dataset.data.data[80, 18, :])
    plt.show()

print(tso.observation.dataset.isBackgroundSubtracted)
assert tso.observation.dataset.isBackgroundSubtracted is True

# sigma clip data
tso.execute("sigma_clip_data")

plt.imshow(np.ma.median(tso.observation.dataset.data[:, :, :], axis=2),
           origin='lower',
           cmap='hot',
           interpolation='nearest',
           aspect='equal')
plt.xlabel("Pixel number spatial direction")
plt.ylabel("Pixel number")
plt.title('Median science data')
plt.show()

plt.imshow(tso.observation.dataset.mask[:, :, :].any(axis=2),
           origin='lower',
           cmap='hot',
           interpolation='nearest',
           aspect='equal')
plt.xlabel("Pixel number spatial direction")
plt.ylabel("Pixel number")
plt.title('Collapsed mask science data')
plt.show()

time0 = tso.observation.dataset.time[80, 18, :].copy()
data0 = tso.observation.dataset.data[80, 18, :].copy()
time0.set_fill_value(np.nan)
data0.set_fill_value(np.nan)
with quantity_support():
    plt.plot(time0.filled(),
             data0.filled())
    plt.show()

print(tso.observation.dataset.isSigmaCliped)
assert tso.observation.dataset.isSigmaCliped is True

# eclipse model
tso.execute("define_eclipse_model")

plt.imshow(tso.model.light_curve_interpolated[0][:, 0, :],
           origin='lower',
           cmap='Reds',
           interpolation='nearest')
plt.colorbar().set_label("Normalised depth")
plt.xlabel("Number of integration")
plt.ylabel("Pixel number dispersion direction")
plt.title('Lightcurve model')
plt.show()

with quantity_support():
    plt.plot(tso.observation.dataset.time.data[80, 18, :],
             tso.model.light_curve_interpolated[0][80, 18, :])
    plt.xlabel("Phase")
    plt.ylabel("Normalised depth")
    plt.show()

# calibration signal
plt.imshow(tso.model.calibration_signal[0][:, 0, :],
           origin='lower',
           cmap='Reds',
           interpolation='nearest')
plt.colorbar().set_label("Normalised depth")
plt.xlabel("Number of integration")
plt.ylabel("Pixel number dispersion direction")
plt.title('Calibration lightcurve model')
plt.show()

print(tso.model.transit_timing)
assert tso.model.transit_timing == [0.4827232723272327, 0.5172767276727672]

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

plt.plot(tso.observation.dataset.time.data.value[80, 18, :],
         tso.cpm.position[80, 18, :])
plt.show()

plt.plot(tso.observation.spectral_trace['wavelength_pixel'].value,
         tso.observation.spectral_trace['positional_pixel'].value -
         (tso.cpm.spectral_trace+tso.cpm.median_position))
axes = plt.gca()
axes.set_ylim([-1, 1])
plt.show()

from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans
from skimage.feature import register_translation
image0 = np.ma.array(tso.observation.dataset.data.data.value[:, 10:30, 0],
                     mask = tso.observation.dataset.mask[:, 10:30, 0])
np.ma.set_fill_value(image0, float("NaN"))
kernel = Gaussian2DKernel(x_stddev=0.3, y_stddev=2.0, theta=-0.1)
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
    image = np.ma.array(tso.observation.dataset.data.data.value[:, 10:30, it],
                     mask = tso.observation.dataset.mask[:, 10:30, it])
    np.ma.set_fill_value(image, float("NaN"))
    cleaned_image = interpolate_replace_nans(image.filled(),
                                         kernel)
    shift, error, diffphase = register_translation(cleaned_image0, cleaned_image, 200)
    shift_store[:,it] = shift

fig = plt.figure(figsize=(8, 5))
ax1 = plt.subplot(1, 2, 1, adjustable='box-forced')
ax2 = plt.subplot(1, 2, 2, sharex=ax1, sharey=ax1, adjustable='box-forced')
ax1.plot(tso.observation.dataset.time.data.value[80, 18, :], shift_store[0,:])
ax1.set_title('Y offset')
ax2.plot(tso.observation.dataset.time.data.value[80, 18, :], shift_store[1,:])
ax2.set_title('X offset')
plt.show()

fig = plt.figure(figsize=(8, 5))
ax1 = plt.subplot(1, 1, 1)
ax1.plot(tso.observation.dataset.time.data.value[80, 18, :],
         tso.cpm.position[80, 18, :])
#ax1.plot(tso.observation.dataset.time.data.value[80, 18, :],-shift_store[0,:] - np.median(-shift_store[0,:]))
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
print(len(tso.cpm.regressor_list))
print(len(tso.cpm.regressor_list[0]))
print(len(tso.cpm.regressor_list[0][50]))
print(len(tso.cpm.regressor_list[0][50][1]))
print(len(tso.cpm.regressor_list[0][50][1][0]))
print(tso.cpm.regressor_list[0][50])

# setup of regression matrix
tso.return_all_design_matrices(center_matrix=True)
print(tso.cpm.design_matrix[0][0][0].shape)
mask = tso.cpm.design_matrix[0][50][0].mask
if mask.shape != ():
    plt.imshow(mask)
    plt.show()
data_masked = tso.cpm.design_matrix[0][50][0]
plt.imshow(data_masked)
plt.show()

B, lambdas = cascade.cpm_model.return_PCR(data_masked.data.T, 36)
plt.plot(B[:, 0:36])
plt.show()

clip_pctl_regressors = float(tso.cascade_parameters.cpm_clip_percentile_regressors)
clip_pctl_time = float(tso.cascade_parameters.cpm_clip_percentile_time)
tso.return_all_design_matrices(clip=True,
                               clip_pctl_time=clip_pctl_time,
                               clip_pctl_regressors=clip_pctl_regressors,
                               center_matrix=False)
print(tso.cpm.design_matrix[0][0][0].shape)
mask = tso.cpm.design_matrix[0][50][0].mask
plt.imshow(mask)
plt.show()
data_masked = tso.cpm.design_matrix[0][50][0]
plt.imshow(data_masked)
plt.show()

B, lambdas = cascade.cpm_model.return_PCR(data_masked.data.T, 36)
plt.plot(B[:, 0:36])
plt.show()

# create calibrated time series and derive planetary signal
tso.execute("calibrate_timeseries")
print(np.ma.max(tso.calibration_results.regularization))
print(np.ma.min(tso.calibration_results.regularization))
plt.imshow(tso.calibration_results.regularization)
plt.show()

# extract planetary signal
tso.execute("extract_spectrum")

extent = [0, 127, 15.5, 7.5]
delta_lam = 15.5-7.5
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
plt.show()

path_old = '/home/bouwman/SST_OBSERVATIONS/projects_HD189733/REDUCED_DATA/'
spec_instr_model_ian = ascii.read(path_old+'results_ian.dat', data_start=1)
wave_ian = (spec_instr_model_ian['lam_micron']*u.micron)
flux_ian = (spec_instr_model_ian['depthA'] *
            u.dimensionless_unscaled).to(u.percent)
error_ian = (spec_instr_model_ian['errorA'] *
             u.dimensionless_unscaled).to(u.percent)
fig, ax = plt.subplots(figsize=(7, 4))
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
ax.plot(wave_ian, flux_ian, color="r", lw=2, alpha=0.9)
ax.errorbar(wave_ian.value, flux_ian.value, yerr=error_ian.value,
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
axes.set_xlim([7.5, 15.5])
axes.set_ylim([0.00, 0.8])
ax.set_ylabel('Fp/Fstar')
ax.set_xlabel('Wavelength')
plt.show()


for i in range(71):
    ks, prob = stats.ks_2samp(flux_ian.value,
                         (tso.exoplanet_spectrum.spectrum.data.data.
                          value[~tso.exoplanet_spectrum.spectrum.mask])[::-1] -
                          (0.035-i*0.001))
    print(-(0.035-i*0.001), ks, prob)


# save planetary signal
tso.execute("save_results")

# plot planetary signal
tso.execute("plot_results")

cal_signal_depth = float(tso.cascade_parameters.cpm_calibration_signal_depth)
mean_eclipse_depth = float(tso.cascade_parameters.observations_median_signal)

plt.plot(tso.calibration_results.calibrated_time_series[90, 17:21, :].T)
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

plt.plot(tso.exoplanet_spectrum.calibration_correction.wavelength,
         tso.exoplanet_spectrum.calibration_correction.data)
#plt.errorbar(tso.exoplanet_spectrum.spectrum.wavelength,
#             tso.exoplanet_spectrum.calibration_correction.data,
#             yerr=tso.exoplanet_spectrum.calibration_correction.uncertainty)
plt.show()

mean_ecld = np.ma.array((mean_eclipse_depth*u.dimensionless_unscaled).to(u.percent))
plt.plot(tso.exoplanet_spectrum.spectrum.wavelength,
         (mean_ecld-tso.exoplanet_spectrum.calibration_correction.data) /
         tso.exoplanet_spectrum.calibration_correction.uncertainty)
plt.plot(tso.exoplanet_spectrum.spectrum.wavelength,
         tso.exoplanet_spectrum.spectrum.data/tso.exoplanet_spectrum.spectrum.uncertainty)
plt.show()
