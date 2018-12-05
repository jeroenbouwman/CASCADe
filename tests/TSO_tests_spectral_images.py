# -*- coding: utf-8 -*-
"""
CASCADe

Test reading Spitzer spectroscopic images

@author: Jeroen Bouwman
"""
import cascade
import os
from matplotlib import pyplot as plt
from astropy.visualization import quantity_support
import numpy as np
from astropy.io import ascii
from astropy import units as u
import pandas as pd
from pandas.plotting import autocorrelation_plot
import seaborn as sns
import warnings
from scipy.io import readsav
from skimage.feature import register_translation
# from scipy.linalg import pinv2
from time import sleep

warnings.simplefilter("ignore")

sns.set_style("white")
# sns.set_context("talk")
sns.set_context("notebook", font_scale=1.5, rc={"lines.linewidth": 2.5})

# create transit spectoscopy object
tso = cascade.TSO.TSOSuite()

# initialization with ini files
path = cascade.initialize.default_initialization_path
tso = cascade.TSO.TSOSuite("cascade_test_cpm2.ini",
                           "cascade_test_object.ini",
                           "cascade_test_data_spectral_images2.ini", path=path)
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
tso.execute("initialize", "cascade_test_cpm2.ini", "cascade_test_object.ini",
            "cascade_test_data_spectral_images2.ini", path=path)
print(tso.cascade_parameters)
print(cascade.initialize.cascade_configuration)
print(tso.cascade_parameters.isInitialized)
assert tso.cascade_parameters == cascade.initialize.cascade_configuration
assert tso.cascade_parameters.isInitialized is True
sleep(1.0)

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

wave0 = tso.observation.dataset.wavelength
wave0_min = np.ma.min(wave0).data.value
wave0_max = np.ma.max(wave0).data.value

med_science_data_image = \
    np.ma.median(tso.observation.dataset.data[:, :, :], axis=2)
med_science_data_image.set_fill_value(np.nan)
plt.imshow(med_science_data_image.filled().value,
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

# create a cleaned version of the spectral data
tso.execute("create_cleaned_dataset")

extracted_spectra = np.ma.sum(tso.cpm.cleaned_data, axis=1)
extracted_spectra.set_fill_value(np.nan)
with quantity_support():
    plt.plot(extracted_spectra.filled())
    plt.title("1D spectra based on cleaned data")
    plt.show()


# determine position of source from data set
tso.execute("determine_source_position")

assert hasattr(tso.cpm, 'spectral_trace') is True
assert hasattr(tso.cpm, 'position') is True

with quantity_support():
    plt.plot(tso.observation.spectral_trace['wavelength_pixel'],
             tso.observation.spectral_trace['positional_pixel'], lw=3)
    plt.title("Shift in spatial diraction of the spectral trace")
    plt.show()
with quantity_support():
    plt.plot(tso.observation.spectral_trace['wavelength_pixel'],
             tso.observation.spectral_trace['wavelength'], lw=3)
    plt.title("Wavelength assignment spectral trace")
    plt.show()

plt.plot(tso.cpm.spectral_trace)
plt.title("Spectral Trace")
plt.show()
plt.plot(tso.cpm.position[80, 18, :])
plt.title("Relative position in slit")
plt.show()

###############################################
# check derived pointing and movement of source
###############################################

roi_mask = tso.observation.instrument_calibration.roi.copy()
kernel = tso.observation.instrument_calibration.convolution_kernel

image_test_roi = np.ma.array(tso.observation.dataset.data.data.value.copy(),
                             mask=tso.observation.dataset.mask.copy())
image_test_roi[roi_mask] = 0.0
image_test_roi.mask[roi_mask] = False
image_test_roi.set_fill_value(np.nan)

cleaned_image_test_roi = tso.cpm.cleaned_data

it0 = np.argmin(np.sum(image_test_roi.mask, axis=(0, 1)))
plt.imshow(image_test_roi[..., it0])
plt.show()
plt.imshow(np.ma.array(cleaned_image_test_roi[..., it0].data.value,
                       mask=cleaned_image_test_roi[..., it0].mask))
plt.show()

_, _, nintegrations = tso.observation.dataset.data.data.value.shape
shift_store = np.zeros((2, nintegrations))
for it in range(nintegrations):
    # subpixel precision
    shift, error, diffphase = \
        register_translation(cleaned_image_test_roi[..., it0],
                             cleaned_image_test_roi[..., it], 200)
    shift_store[:, it] = shift

fig = plt.figure(figsize=(8, 5))
ax1 = plt.subplot(1, 2, 1, adjustable='box')
ax2 = plt.subplot(1, 2, 2, sharex=ax1, sharey=ax1, adjustable='box')
ax1.plot(tso.observation.dataset.time.data.value[80, 18, :], shift_store[0, :])
ax1.set_title('Y offset')
ax2.plot(tso.observation.dataset.time.data.value[80, 18, :], shift_store[1, :])
ax2.set_title('X offset')
plt.show()

fig = plt.figure(figsize=(10, 6))
ax1 = plt.subplot(1, 1, 1)
ax1.plot(tso.observation.dataset.time.data.value[80, 18, :],
         tso.cpm.position[80, 18, :], label='CoL, lw=2')
ax1.plot(tso.observation.dataset.time.data.value[80, 18, :],
         shift_store[0, :]-np.median(shift_store[0, :]), label='CC disp')
ax1.plot(tso.observation.dataset.time.data.value[80, 18, :],
         -shift_store[1, :]-np.median(-shift_store[1, :]), label='CC cross',
         lw=2)
ax1.legend(loc='best')
ax1.set_ylim([-0.05, 0.05])
plt.show()

# read file containing derived pointing using peakup array
path_to_save_file = \
    ('/home/bouwman/SST_OBSERVATIONS/projects_HD189733/REDUCED_DATA/'
     'HD_189733_combined_SL_B/S16.1.0/SL1/')
filename_save_file = 'pointing_offsets_sl1_B.save'
movement_AOR = readsav(os.path.join(path_to_save_file, filename_save_file))

xmovement_slit = movement_AOR['x_offset_spectrum_fit_sl1_b']
xmovement_peakup = movement_AOR['x_offset_image_fit_sl1_b']
ymovement_peakup = movement_AOR['offset_image_fit_sl1_b']
phase = movement_AOR['phase_spectrum_fit_sl1_b']

fig = plt.figure(figsize=(10, 8))
ax1 = plt.subplot(1, 1, 1)
ax1.plot(tso.observation.dataset.time.data.value[80, 18, :],
         tso.cpm.position[80, 18, :], label='COL')
ax1.plot(tso.observation.dataset.time.data.value[80, 18, :],
         -shift_store[1, :]-np.median(-shift_store[1, :]), label='CC')
ax1.plot(phase, xmovement_slit-np.median(xmovement_slit), label='Slit')
ax1.plot(phase, xmovement_peakup-np.median(xmovement_peakup), label='PU')
ax1.set_ylim([-0.05, 0.05])
ax1.legend(loc='best')
ax1.set_title("Comparison cross-dispersion direction")
plt.show()

fig = plt.figure(figsize=(10, 8))
ax1 = plt.subplot(1, 1, 1)
# ax1.plot(tso.observation.dataset.time.data.value[80, 18, :],
#         tso.cpm.position[80, 18, :], label='COL')
ax1.plot(tso.observation.dataset.time.data.value[80, 18, :],
         (-shift_store[0, :]-np.median(-shift_store[0, :])), label='CC')
ax1.plot(phase, ymovement_peakup-np.median(ymovement_peakup), label='PU')
ax1.set_ylim([-0.2, 0.2])
ax1.legend(loc='best')
ax1.set_title("Comparison dispersion direction")
plt.show()

#######################################################

# set the extraction area
tso.execute("set_extraction_mask")

print(hasattr(tso.cpm, 'extraction_mask'))
assert hasattr(tso.cpm, 'extraction_mask') is True
print(tso.cpm.extraction_mask[0].shape)
assert (tso.cpm.extraction_mask[0].shape[0] ==
        tso.observation.dataset.data.shape[0])

plt.imshow(tso.cpm.extraction_mask[0],
           origin='lower',
           cmap='Reds',
           interpolation='nearest')
plt.ylabel("Data number in wavelength direction")
plt.xlabel("Data number in spatial direction")
plt.title("Extraction Mask")
plt.show()

# optimally extract spectrum of target star
tso.execute("optimal_extraction")

assert hasattr(tso.cpm, 'extraction_profile') is True

extracted_spectra = tso.observation.dataset_optimal_extracted.data
time = tso.observation.dataset_optimal_extracted.time
wavelength = tso.observation.dataset_optimal_extracted.wavelength
cleaned_data = tso.cpm.cleaned_data
cleaned_data.set_fill_value(np.NaN)

mean_wavelength = np.ma.mean(wavelength, axis=1)
mean_wavelength.set_fill_value(np.NaN)
mean_time = np.ma.mean(time, axis=0)
mean_time.set_fill_value(np.NaN)
mean_spectra2 = np.ma.mean(extracted_spectra, axis=1)
mean_spectra2.set_fill_value(np.NaN)

with quantity_support():
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(mean_time.filled(),
            np.ma.mean(np.ma.sum(cleaned_data, axis=1), axis=0),
            color='gray', lw=3, label='Tapered')
    ax.plot(mean_time.filled(),
            np.ma.mean(extracted_spectra, axis=0).filled(), lw=3,
            alpha=0.6, color='red', label='Optimal')
    ax.legend(loc='best')
    ax.set_title('Comparison optimal and tapered extraction')
    ax.set_xlabel('Phase')
    plt.show()
with quantity_support():
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(mean_wavelength.filled(),
            np.ma.mean(np.ma.sum(cleaned_data, axis=1).filled(), axis=1),
            color='gray', lw=3, label='Tapered')
    ax.plot(mean_wavelength.filled(),
            mean_spectra2.filled(), lw=3,
            alpha=0.6, color='red', label='Optimal')
    ax.legend(loc='best')
    ax.set_title('Comparison optimal and tapered extraction')
    plt.show()

# setup regressors
tso.execute("select_regressors")

print(hasattr(tso.cpm, 'regressor_list'))
assert hasattr(tso.cpm, 'regressor_list') is True

# setup of regression matrix without clipping
tso.return_all_design_matrices()

# 2D case:
#   regressor list index:
#       first: [# nod]
#       second: [# of valid pixels within extraction mask]
#       third: [0=pixel coord; 1=list of regressors]
#       forth: [0=coordinate wave direction; 1=coordinate spatial direction]
#   design_matrix:
#       first: [# nod]
#       second: [# of valid pixels within extraction mask]
#       third: [0]


i_pixel_in_extraction_mask = 100
reg_list = \
    tso.cpm.regressor_list[0][i_pixel_in_extraction_mask][1][0]
reg_matrix = tso.cpm.design_matrix[0][i_pixel_in_extraction_mask][0]

plt.imshow(reg_matrix.mask,
           origin='lower',
           cmap='hot',
           interpolation='nearest',
           aspect='auto')
plt.xlabel("Integration Number")
plt.ylabel("Regressor Number")
plt.title('Mask of non-clipped regressor matrix')
plt.show()

data_masked = reg_matrix
data_masked.set_fill_value(np.nan)
plt.imshow(data_masked.filled().value,
           origin='lower',
           cmap='hot',
           interpolation='nearest',
           aspect='auto')
plt.colorbar().set_label("Intensity ({})".format(image0.data.unit))
plt.xlabel("Integration Number")
plt.ylabel("Regressor Number")
plt.title('Not clipped regressor matrix')
plt.show()

# setup of regression matrix with clipping
clip_pctl_time = float(tso.cascade_parameters.cpm_clip_percentile_time)
clip_pctl_regressors = float(tso.cascade_parameters.
                             cpm_clip_percentile_regressors)
tso.return_all_design_matrices(clip=True, clip_pctl_time=clip_pctl_time,
                               clip_pctl_regressors=clip_pctl_regressors)

reg_matrix = tso.cpm.design_matrix[0][i_pixel_in_extraction_mask][0]

plt.imshow(reg_matrix.mask,
           origin='lower',
           cmap='hot',
           interpolation='nearest',
           aspect='auto')
plt.xlabel("Integration Number")
plt.ylabel("Regressor Number")
plt.title('Mask of clipped regressor matrix')
plt.show()

data_masked = reg_matrix
data_masked.set_fill_value(np.nan)
plt.imshow(data_masked.filled().value,
           origin='lower',
           cmap='hot',
           interpolation='nearest',
           aspect='auto')
plt.colorbar().set_label("Intensity ({})".format(reg_matrix.data.unit))
plt.xlabel("Integration Number")
plt.ylabel("Regressor Number")
plt.title('Clipped regressor matrix')
plt.show()

plt.plot(data_masked.filled().value[30:40, :].T)
plt.show()

# eclipse model
tso.execute("define_eclipse_model")

plt.imshow(tso.model.light_curve_interpolated[0][:, 0, :],
           origin='lower',
           cmap='Reds',
           interpolation='nearest',
           aspect='auto',
           extent=[0, image0.shape[1], wave0_min, wave0_max])
plt.colorbar().set_label("Normalised depth")
plt.xlabel("Number of integration")
plt.ylabel("Wavelength")
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
sleep(1.0)

# create calibrated time series and derive planetary signal
tso.execute("calibrate_timeseries")

print(np.ma.max(tso.calibration_results.regularization))
print(np.ma.min(tso.calibration_results.regularization))
assert (np.ma.max(tso.calibration_results.regularization) <
        float(tso.cascade_parameters.cpm_lam1))
assert (np.ma.min(tso.calibration_results.regularization) >
        float(tso.cascade_parameters.cpm_lam0))

fig, ax = plt.subplots(figsize=(7, 5))
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
im = ax.imshow(tso.calibration_results.regularization,
               origin='lower',
               cmap='hot',
               interpolation='nearest',
               aspect='auto')
fig.colorbar(im, ax=ax).set_label("Regularization strength")
ax.set_ylabel("Wavelength data point number")
ax.set_xlabel("Spatial data point number")
ax.set_title('Regularization parameter')
plt.show()

with quantity_support():
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(tso.observation.dataset.time.data[80, 18, :],
            tso.observation.dataset.data.data[80, 18, :], label='Data')
    ax.plot(tso.observation.dataset.time.data[80, 18, :],
            tso.calibration_results.model_time_series.data[80, 18, :],
            label='Model')
    ax.legend(loc='best')
    ax.set_title('Comparison Data and Regression Model')
    ax.set_xlabel('Phase')
    plt.show()


# extract planetary signal
tso.execute("extract_spectrum")

extent = [0, tso.exoplanet_spectrum.weighted_image.shape[1]-1,
          wave0_min, wave0_max]
delta_lam = wave0_max - wave0_min
fig, ax = plt.subplots(figsize=(10, 10))
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

# correct the extracted planetary signal for non uniform subtraction
# of averige signal
tso.execute("correct_extracted_spectrum")

wavelength_corrected_spectrum = \
    tso.exoplanet_spectrum.corrected_spectrum.wavelength
corrected_spectrum = tso.exoplanet_spectrum.corrected_spectrum.data
error_corrected_spectrum = \
    tso.exoplanet_spectrum.corrected_spectrum.uncertainty

path_old = '/home/bouwman/SST_OBSERVATIONS/projects_HD189733/REDUCED_DATA/'
spec_instr_model_ian = ascii.read(path_old+'results_ian.dat', data_start=1)
wave_ian = (spec_instr_model_ian['lam_micron']*u.micron)
flux_ian = (spec_instr_model_ian['depthB'] *
            u.dimensionless_unscaled).to(u.percent)
error_ian = (spec_instr_model_ian['errorB'] *
             u.dimensionless_unscaled).to(u.percent)
fig, ax = plt.subplots(figsize=(13, 9))
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
ax.plot(wavelength_corrected_spectrum, corrected_spectrum,
        lw=3, zorder=6, color='green')
ax.errorbar(wavelength_corrected_spectrum.data.value,
            corrected_spectrum.data.value,
            yerr=error_corrected_spectrum.data.value,
            fmt=".k", color='green', lw=3, alpha=0.5, ecolor='green',
            markerfacecolor='green',
            markeredgecolor='green', fillstyle='full', markersize=10,
            zorder=6)
axes = plt.gca()
axes.set_xlim([7.5, 15.5])
axes.set_ylim([0.00, 0.8])
ax.set_ylabel('Fp/Fstar')
ax.set_xlabel('Wavelength')
plt.show()

image0 = tso.exoplanet_spectrum.weighted_residual.copy()
image_unit = image0.data.unit
image0 = np.ma.array(image0.data,
                     mask=image0.mask)
image0.set_fill_value(np.nan)
image0 = image0.filled().value

wave0 = tso.observation.dataset.wavelength
wave0_min = np.ma.min(wave0).data.value
wave0_max = np.ma.max(wave0).data.value
plt.imshow(image0,
           origin='lower',
           cmap='hot',
           interpolation='nearest',
           aspect='auto',
           extent=[0, image0.shape[1], wave0_min, wave0_max])
plt.colorbar().set_label("Residual ({})".format(image_unit))
plt.xlabel("Integration Number")
plt.ylabel("Wavelength ")
plt.title('Residual Image')
plt.show()

df = pd.DataFrame(image0.T)
for col_id in df.dropna(axis=1, how='all').columns:
    sns.distplot(df.dropna(axis=1, how='all')[col_id].dropna())


wave_spectrum = tso.exoplanet_spectrum.spectrum.wavelength.copy()
wave_spectrum.set_fill_value(np.nan)
wave_spectrum = wave_spectrum.filled().value
fig, ax = plt.subplots(figsize=(7, 5))
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
ax.plot(wave_spectrum, tso.exoplanet_spectrum.weighted_aic,
        label='AIC')
ax.legend(loc='best')
ax.set_title('AIC model fit per wavelength')
ax.set_xlabel('Wavelength')
ax.set_xlabel('AIC')
plt.show()

df = pd.DataFrame(image0.T)
df = df.replace(np.nan, 0)
fig, ax = plt.subplots(figsize=(7, 5))
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
ax = autocorrelation_plot(df[10], lw=3, ax=ax, label="10")
ax = autocorrelation_plot(df[50], lw=3, ax=ax, label="50")
ax = autocorrelation_plot(df[90], lw=3, ax=ax, label="90")
ax = autocorrelation_plot(df, lw=3, ax=ax, label="All")
ax.legend(loc='best')
ax.set_title('ACF in time direction')
plt.show()

df = pd.DataFrame(image0)
df = df.replace(np.nan, 0)
fig, ax = plt.subplots(figsize=(7, 5))
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
ax = autocorrelation_plot(df[50], lw=3, ax=ax, label="50")
ax = autocorrelation_plot(df[100], lw=3, ax=ax, label="100")
ax = autocorrelation_plot(df[150], lw=3, ax=ax, label="150")
ax = autocorrelation_plot(df[250], lw=3, ax=ax, label="250")
ax = autocorrelation_plot(df, lw=3, ax=ax, label="All")
ax.legend(loc='best')
ax.set_title('ACF in wavelength direction')
plt.show()


fig, ax = plt.subplots(figsize=(7, 5))
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
s = pd.Series(np.random.normal(0.0, 1.0, 280))
ax = autocorrelation_plot(s, lw=3, ax=ax, label="run1")
s = pd.Series(np.random.normal(0.0, 1.0, 280))
ax = autocorrelation_plot(s, lw=3, ax=ax, label="run2")
s = pd.Series(np.random.normal(0.0, 1.0, 280))
ax = autocorrelation_plot(s, lw=3, ax=ax, label="run3")
s = pd.Series(np.random.normal(0.0, 1.0, 280))
ax = autocorrelation_plot(s, lw=3, ax=ax, label="run4")
s = pd.DataFrame(np.random.normal(0.0, 1.0, 280*130).reshape(280, 130))
ax = autocorrelation_plot(s, lw=3, ax=ax, label="All")
ax.legend(loc='best')
ax.set_title('ACF in time direction for random distributions')
plt.show()

fig, ax = plt.subplots(figsize=(7, 5))
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
s = pd.Series(np.random.normal(0.0, 1.0, 130))
ax = autocorrelation_plot(s, lw=3, ax=ax, label="run1")
s = pd.Series(np.random.normal(0.0, 1.0, 130))
ax = autocorrelation_plot(s, lw=3, ax=ax, label="run2")
s = pd.Series(np.random.normal(0.0, 1.0, 130))
ax = autocorrelation_plot(s, lw=3, ax=ax, label="run3")
s = pd.Series(np.random.normal(0.0, 1.0, 130))
ax = autocorrelation_plot(s, lw=3, ax=ax, label="run4")
s = pd.DataFrame(np.random.normal(0.0, 1.0, 280*130).reshape(130, 280))
ax = autocorrelation_plot(s, lw=3, ax=ax, label="All")
ax.legend(loc='best')
ax.set_title('ACF in wavelength direction for random distributions')
plt.show()

# save planetary signal
tso.execute("save_results")

# plot planetary signal
tso.execute("plot_results")
