# -*- coding: utf-8 -*-

import numpy as np
from matplotlib import pyplot as plt
from astropy.visualization import quantity_support
from astropy.io import ascii
from astropy import units as u
import pandas as pd
from pandas.plotting import autocorrelation_plot
import cascade
from cascade.utilities import spectres
import seaborn as sns

sns.set_style("white")

# create transit spectoscopy object
tso = cascade.TSO.TSOSuite()

# reset
tso.execute("reset")

# initialization with ini files
path = cascade.initialize.default_initialization_path
tso = cascade.TSO.TSOSuite("initialize", "cascade_test_HST_cpm2.ini",
                           "cascade_test_HST_object2.ini",
                           "cascade_test_HST_data_spectral_images2.ini",
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

assert tso.observation.dataset.time.data.shape == (256, 256, 484)
assert type(tso.observation.dataset.time.data) == u.quantity.Quantity
assert tso.observation.dataset.time.data.unit == u.dimensionless_unscaled

# subtract background
tso.execute("subtract_background")

# sigma clip data
tso.execute("sigma_clip_data")

wave0 = tso.observation.dataset.wavelength
wave0_min = np.ma.min(wave0).data.value
wave0_max = np.ma.max(wave0).data.value

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

time0 = tso.observation.dataset.time[80, 85, :].copy()
data0 = tso.observation.dataset.data[80, 85, :].copy()
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
    plt.plot(tso.observation.dataset.time.data[80+64, 85+64, :],
             tso.model.light_curve_interpolated[0][80+64, 85+64, :])
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
plt.plot(tso.cpm.position[80+64, 85+64, :])
plt.title("Relative position in slit")
plt.show()


###################################
from skimage.feature import register_translation

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


def pad_with(vector, pad_width, iaxis, kwargs):
    pad_value = kwargs.get('padder', 10)
    vector[:pad_width[0]] = pad_value
    vector[-pad_width[1]:] = pad_value
    return vector
# cleaned_image_test0 = np.pad(cleaned_image_test0, 2, pad_with, padder=0)


_, _, nintegrations = tso.observation.dataset.data.data.value.shape
shift_store = np.zeros((2, nintegrations))
for it in range(nintegrations):
    # subpixel precision
    shift, error, diffphase = \
        register_translation(cleaned_image_test_roi[..., it0],
                             cleaned_image_test_roi[..., it], 200)
    shift_store[:, it] = shift

fig = plt.figure(figsize=(9, 6))
ax1 = plt.subplot(1, 2, 1, adjustable='box')
ax2 = plt.subplot(1, 2, 2, sharex=ax1, sharey=ax1, adjustable='box')
ax1.plot(tso.observation.dataset.time.data.value[80+64, 85+64, :],
         shift_store[0, :])
ax1.set_title('Y offset')
ax2.plot(tso.observation.dataset.time.data.value[80+64, 85+64, :],
         shift_store[1, :])
ax2.set_title('X offset')
plt.show()

fig = plt.figure(figsize=(9, 6))
ax1 = plt.subplot(1, 1, 1)
ax1.plot(tso.observation.dataset.time.data.value[80+64, 85+64, :],
         tso.cpm.position[80+64, 85+64, :], label='CoL', lw=2)
ax1.plot(tso.observation.dataset.time.data.value[80, 18, :],
         shift_store[0, :]-np.median(shift_store[0, :]), label='CC disp')
ax1.plot(tso.observation.dataset.time.data.value[80+64, 18, :],
         -shift_store[1, :]-np.median(-shift_store[1, :]), label='CC cross',
         lw=2)
ax1.legend(loc='best')
ax1.set_ylim([-0.08, 0.08])
plt.show()

fig = plt.figure(figsize=(9, 6))
ax1 = plt.subplot(1, 1, 1)
ax1.plot(tso.observation.dataset.time.data.value[80+64, 85+64, :],
         tso.cpm.position[80+64, 85+64, :], label='COL')
ax1.plot(tso.observation.dataset.time.data.value[80+64, 85+64, :],
         -shift_store[1, :]-np.median(-shift_store[1, :]), label='CC')
ax1.set_ylim([-0.08, 0.08])
ax1.legend(loc='best')
ax1.set_title("Comparison cross-dispersion direction")
plt.show()

fig = plt.figure(figsize=(9, 6))
ax1 = plt.subplot(1, 1, 1)
# ax1.plot(tso.observation.dataset.time.data.value[80, 18, :],
#         tso.cpm.position[80, 18, :], label='COL')
ax1.plot(tso.observation.dataset.time.data.value[80+64, 85+64, :],
         (-shift_store[0, :]-np.median(-shift_store[0, :])), label='CC')
ax1.set_ylim([-0.1, 0.1])
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


i_pixel_in_extraction_mask = 150
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
plt.colorbar().set_label("Intensity ({})".format(image0.data.unit))
plt.xlabel("Integration Number")
plt.ylabel("Regressor Number")
plt.title('Clipped regressor matrix')
plt.show()

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
    ax.plot(tso.observation.dataset.time.data[80+64, 151, :],
            tso.observation.dataset.data.data[80+64, 151, :], label='Data')
    ax.plot(tso.observation.dataset.time.data[80+64, 151, :],
            tso.calibration_results.model_time_series.data[80+64, 151, :],
            label='Model')
    ax.legend(loc='best')
    ax.set_title('Comparison Data and Regression Model')
    ax.set_xlabel('Phase')
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
spec_instr_model_mandell = ascii.read(path_old+'WASP12b_mandell.txt',
                                      data_start=0)
wave_mandell = (spec_instr_model_mandell['col1']*u.micron)
flux_mandell = (spec_instr_model_mandell['col2']*u.percent)
error_mandell = (spec_instr_model_mandell['col3']*u.percent)

path_old = '/home/bouwman//'
spec_instr_model_tsiaras = ascii.read(path_old+'wasp12b_tsiaras.txt',
                                      data_start=0)
wave_tsiaras = (spec_instr_model_tsiaras['col1']*u.micron)
flux_tsiaras = (spec_instr_model_tsiaras['col2']*100*u.percent)

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

W = (tso.exoplanet_spectrum.weighted_normed_parameters[:-4, :] /
     np.ma.sum(tso.exoplanet_spectrum.weighted_normed_parameters[:-1, :],
               axis=0)).T
K = W - np.identity(W.shape[0])
from scipy.linalg import pinv2
K.set_fill_value(0.0)
weighted_signal = tso.exoplanet_spectrum.weighted_signal.copy()
weighted_signal.set_fill_value(0.0)

bla = np.dot(pinv2(K.filled(), rcond=0.8),
             -weighted_signal.filled())


from scipy.linalg import lstsq
res = lstsq(K.filled(), -weighted_signal.filled(), cond=0.8)
bla=res[0]


reg_par = {'lam0': 1.e-7, 'lam1': 1.e-2, 'nlam': 100}
lam_reg0 = reg_par["lam0"]  # lowest value of regularization parameter
lam_reg1 = reg_par["lam1"]   # highest
ngrid_lam = reg_par["nlam"]  # number of points in grid
# array to hold values of regularization parameter grid
delta_lam = np.abs(np.log10(lam_reg1) - np.log10(lam_reg0)) / (ngrid_lam-1)
lam_reg_array = 10**(np.log10(lam_reg0) +
                     np.linspace(0, ngrid_lam-1, ngrid_lam)*delta_lam)

from sklearn import linear_model
RCV = linear_model.RidgeCV(fit_intercept=False, alphas=lam_reg_array)
#RCV = linear_model.RidgeCV(fit_intercept=False)
RCV.fit(K.filled(), -weighted_signal.filled())
bla = RCV.coef_

bla = bla-np.ma.median(bla)
planet_radius = \
    (u.Quantity(tso.cascade_parameters.object_radius).to(u.m) /
     u.Quantity(tso.cascade_parameters.
                object_radius_host_star).to(u.m))
planet_radius = planet_radius.decompose().value
bla = (bla * (1.0 - planet_radius**2) +
       planet_radius**2)

fig, ax = plt.subplots(figsize=(7, 4))
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
ax.plot(rebinned_wave, rebinned_spec, color="black", lw=4, alpha=0.9)
ax.errorbar(rebinned_wave, rebinned_spec, yerr=rebinned_error,
            fmt=".k", color="black", lw=4,
            alpha=0.9, ecolor="black",
            markeredgecolor='black', fillstyle='full', markersize=10,
            markerfacecolor='black', zorder=3)
ax.plot(wave_mandell, flux_mandell, color="r", lw=3, alpha=0.9)
ax.errorbar(wave_mandell.value, flux_mandell.value, yerr=error_mandell.value,
            fmt=".k", color="r", lw=2,
            alpha=0.9, ecolor="r",
            markeredgecolor='r', fillstyle='full', markersize=10,
            markerfacecolor='r', zorder=3)
ax.plot(wave_tsiaras, flux_tsiaras, color="magenta", lw=3, alpha=0.9)
ax.plot(tso.exoplanet_spectrum.spectrum.wavelength,
        tso.exoplanet_spectrum.spectrum.data, lw=3, alpha=0.5, color='blue')
ax.errorbar(tso.exoplanet_spectrum.spectrum.wavelength.data.value,
            tso.exoplanet_spectrum.spectrum.data.data.value,
            yerr=tso.exoplanet_spectrum.spectrum.uncertainty.data.value,
            fmt=".k", color='blue', lw=3, alpha=0.5, ecolor='blue',
            markerfacecolor='blue',
            markeredgecolor='blue', fillstyle='full', markersize=10,
            zorder=4)
ax.plot(tso.exoplanet_spectrum.spectrum.wavelength, bla*100, lw=4, zorder=6)
axes = plt.gca()
axes.set_xlim([1.1, 1.7])
axes.set_ylim([1.2, 1.5])
ax.set_ylabel('Transit depth')
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
ax = autocorrelation_plot(df[10+64], lw=3, ax=ax, label="10")
ax = autocorrelation_plot(df[50+64], lw=3, ax=ax, label="50")
ax = autocorrelation_plot(df[90+64], lw=3, ax=ax, label="90")
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
