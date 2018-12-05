# -*- coding: utf-8 -*-

import cascade
from matplotlib import pyplot as plt
from astropy.visualization import quantity_support
import numpy as np
from astropy.io import ascii
from astropy import units as u
import pandas as pd
from pandas.plotting import autocorrelation_plot
import seaborn as sns
from cascade.utilities import spectres
# from astropy.convolution import Box1DKernel, convolve
# from astropy.stats import sigma_clip
from sklearn.preprocessing import RobustScaler
# from scipy.linalg import pinv2
from time import sleep
import warnings

warnings.simplefilter("ignore")

sns.set_style("white")

# create transit spectoscopy object
tso = cascade.TSO.TSOSuite()

# reset
tso.execute("reset")

# initialization with ini files
path = cascade.initialize.default_initialization_path
tso = cascade.TSO.TSOSuite("initialize",
                           "cascade_test_HST_cpm_spectra2_coe.ini",
                           "cascade_test_HST_object2.ini",
                           "cascade_test_HST_data_spectra2_coe.ini", path=path)

print(tso.cascade_parameters)
print(cascade.initialize.cascade_configuration)
print(tso.cascade_parameters.isInitialized)
assert tso.cascade_parameters == cascade.initialize.cascade_configuration
assert tso.cascade_parameters.isInitialized is True
sleep(1.0)

# load data
tso.execute("load_data")

image0 = tso.observation.dataset.data[:, :].copy()
image0.set_fill_value(np.nan)
wave0 = tso.observation.dataset._wavelength
wave0_min = np.min(wave0)
wave0_max = np.max(wave0)
plt.imshow(image0.filled().value,
           origin='lower',
           cmap='hot',
           interpolation='nearest',
           aspect='auto',
           extent=[0, image0.shape[1], wave0_min, wave0_max])
plt.colorbar().set_label("Intensity ({})".format(image0.data.unit))
plt.xlabel("Integration Number")
plt.ylabel("Wavelength")
plt.title('Science data')
plt.show()


time = tso.observation.dataset.time[80, :].copy()
time.set_fill_value(np.nan)
flux = tso.observation.dataset.data[80, :].copy()
flux.set_fill_value(np.nan)
with quantity_support():
    plt.plot(time.filled(), flux.filled())
    plt.show()
wavelength = tso.observation.dataset.wavelength[:, 0].copy()
wavelength.set_fill_value(np.nan)
flux = tso.observation.dataset.data[:, 0].copy()
flux.set_fill_value(np.nan)
with quantity_support():
    plt.plot(wavelength.filled(), flux.filled())
    plt.show()

print(tso.observation.dataset.time.data.shape)
print(type(tso.observation.dataset.time.data))
print(tso.observation.dataset.time.data.unit)
assert tso.observation.dataset.time.data.shape == (256, 484)
assert type(tso.observation.dataset.time.data) == u.quantity.Quantity
assert tso.observation.dataset.time.data.unit == u.dimensionless_unscaled

# subtract background
tso.execute("subtract_background")

with quantity_support():
    plt.plot(tso.observation.dataset.time[80, :],
             tso.observation.dataset.data[80, :])
    plt.show()

print(hasattr(tso.observation.dataset, 'isBackgroundSubtracted'))
assert hasattr(tso.observation.dataset, 'isBackgroundSubtracted') is False


# sigma clip data
tso.execute("sigma_clip_data")

image0 = tso.observation.dataset.data[:, :].copy()
image0.set_fill_value(np.nan)
wave0 = tso.observation.dataset._wavelength
wave0_min = np.min(wave0)
wave0_max = np.max(wave0)
plt.imshow(image0.filled().value,
           origin='lower',
           cmap='hot',
           interpolation='nearest',
           aspect='auto',
           extent=[0, image0.shape[1], wave0_min, wave0_max])
plt.colorbar().set_label("Intensity ({})".format(image0.data.unit))
plt.xlabel("Integration Number")
plt.ylabel("Wavelength")
plt.title('Science data')
plt.show()

plt.imshow(tso.observation.dataset.mask[:, :],
           origin='lower',
           cmap='hot',
           interpolation='nearest',
           aspect='auto',
           extent=[0, image0.shape[1], wave0_min, wave0_max])
plt.xlabel("Integration Number")
plt.ylabel("Wavelength")
plt.title('Data mask')
plt.show()

time0 = tso.observation.dataset.time[80, :].copy()
data0 = tso.observation.dataset.data[80, :].copy()
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

###################################
roi_mask = tso.observation.instrument_calibration.roi.copy()
kernel = tso.observation.instrument_calibration.convolution_kernel

image_test_roi = np.ma.array(tso.observation.dataset.data.data.value.copy(),
                             mask=tso.observation.dataset.mask.copy())
image_test_roi[roi_mask] = 0.0
image_test_roi.mask[roi_mask] = False
image_test_roi.set_fill_value(np.nan)

cleaned_image_test_roi = tso.cpm.cleaned_data
cleaned_image_test_roi.set_fill_value(np.nan)

plt.imshow(image_test_roi,
           origin='lower',
           cmap='hot',
           interpolation='nearest',
           aspect='auto')
plt.show()
plt.imshow(cleaned_image_test_roi.filled().value,
           origin='lower',
           cmap='hot',
           interpolation='nearest',
           aspect='auto')
plt.show()

RS = RobustScaler(with_scaling=True)
RS.fit(cleaned_image_test_roi.T)
X_scaled_masked = RS.transform(cleaned_image_test_roi.T)
bla = np.ma.array(X_scaled_masked, mask=cleaned_image_test_roi.T.mask)
plt.plot(bla)
plt.plot(np.ma.median(bla, axis=1))
plt.ylim([-4, 4])
plt.show()

RS = RobustScaler(with_scaling=True)
RS.fit(image0.T)
X_scaled_masked = RS.transform(image0.T)
bla = np.ma.array(X_scaled_masked, mask=image0.T.mask)
plt.plot(bla)
plt.plot(np.ma.median(bla, axis=1))
plt.ylim([-4, 4])
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
with quantity_support():
    plt.plot(tso.observation.dataset.time[80, :],
             tso.observation.dataset.position[80, :])
    plt.title("Source position (cross dispersion) in the slit")
    plt.show()

plt.plot(tso.cpm.spectral_trace)
plt.title("Spectral Trace")
plt.show()
plt.plot(tso.cpm.position[80, :])
plt.title("Relative position in slit")
plt.show()

# set the extraction area
tso.execute("set_extraction_mask")

print(hasattr(tso.cpm, 'extraction_mask'))
assert hasattr(tso.cpm, 'extraction_mask') is True
print(tso.cpm.extraction_mask[0].shape)
assert (tso.cpm.extraction_mask[0].shape[0] ==
        tso.observation.dataset.data.shape[0])

plt.plot(tso.cpm.extraction_mask[0], lw=3)
plt.xlabel("Data number in wavelength direction")
plt.ylabel("Data number in spatial direction")
plt.title("Extraction Mask")
plt.show()

# setup regressors
tso.execute("select_regressors")

print(hasattr(tso.cpm, 'regressor_list'))
assert hasattr(tso.cpm, 'regressor_list') is True

# setup of regression matrix without clipping
tso.return_all_design_matrices()

# 1D case:
#   regressor list index:
#       first: [# nod]
#       second: [# of valid pixels within extraction mask]
#       third: [0=pixel coord; 1=list of regressors]
#       forth: [0=coordinate wave direction; 1=coordinate spatial direction]
#   design_matrix:
#       first: [# nod]
#       second: [# of valid pixels within extraction mask]
#       third: [0]

i_data_point_in_wavelength_direction = 50
reg_list = \
    tso.cpm.regressor_list[0][i_data_point_in_wavelength_direction][1][0]
reg_matrix = tso.cpm.design_matrix[0][i_data_point_in_wavelength_direction][0]

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

reg_matrix = tso.cpm.design_matrix[0][i_data_point_in_wavelength_direction][0]

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

# eclipse model
tso.execute("define_eclipse_model")

plt.imshow(tso.model.light_curve_interpolated[0][:, :],
           origin='lower',
           cmap='hot',
           interpolation='nearest',
           aspect='auto',
           extent=[0, image0.shape[1], wave0_min, wave0_max])
plt.xlabel("Integration Number")
plt.ylabel("Wavelength")
plt.title('Lightcurve model')
plt.show()

with quantity_support():
    plt.plot(tso.observation.dataset.time.data[80, :],
             tso.model.light_curve_interpolated[0][80, :])
    plt.show()

RS = RobustScaler(with_scaling=True)
RS.fit(tso.model.light_curve_interpolated[0][80, :].reshape(-1, 1))
X_scaled = \
  RS.transform(tso.model.light_curve_interpolated[0][80, :].reshape(-1, 1))
mean_data = np.ma.mean(tso.cpm.cleaned_data, axis=0)
mean_data.set_fill_value(np.NaN)
RS = RobustScaler(with_scaling=True)
RS.fit(mean_data.filled().reshape(-1, 1))
X2_scaled = \
  RS.transform(mean_data.filled().reshape(-1, 1))

plt.plot(tso.observation.dataset.time.data[80, :], X_scaled.squeeze()*0.88,
         lw=3, alpha=0.6, color='blue', label='Model')
plt.plot(tso.observation.dataset.time.data[80, :], X2_scaled.squeeze(), lw=3,
         alpha=0.6, color='red', label='Data')
plt.xlabel("Phase")
plt.ylabel("Normalised depth")
plt.title('Comparison Model and optimal mean signal')
plt.show()

print(tso.model.transit_timing)
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
ax.plot(tso.calibration_results.regularization, lw=3)
ax.set_xlabel("Wavelength data point number")
ax.set_ylabel("Regularization strength")
ax.set_title('Regularization parameter')
plt.show()

with quantity_support():
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(tso.observation.dataset.time.data[80, :],
            tso.observation.dataset.data.data[80, :], label='Data')
    ax.plot(tso.observation.dataset.time.data[80, :],
            tso.calibration_results.model_time_series.data[80, 0, :],
            label='Model')
    ax.legend(loc='best')
    ax.set_title('Comparison Data and Regression Model')
    ax.set_xlabel('Phase')
    plt.show()


# extract planetary signal
tso.execute("extract_spectrum")

extent = [0, 127, 1.8, 1.05]
delta_lam = 1.8-1.05
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

# correct the extracted planetary signal for non uniform subtraction
# of averige signal
tso.execute("correct_extracted_spectrum")

##############################
# Read resuls from literature
##############################
path_old = '/home/bouwman//'

corrections_wasp12_stevenson = ascii.read(path_old+'WASP12b_corrections.txt',
                                      data_start=1)
wave_correction = (corrections_wasp12_stevenson['col1']*u.micron +
                   corrections_wasp12_stevenson['col2']*u.micron)/2.0
corection_signal = (1.0 +
                    corrections_wasp12_stevenson['col5']+
                    corrections_wasp12_stevenson['col7'])
corrected_signal_stevenson2014 = corrections_wasp12_stevenson['col9']*1.02
error_corrected_signal_stevenson2014 = corrections_wasp12_stevenson['col10']*1.02

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
error_tsiaras = (spec_instr_model_tsiaras['col3']*100*u.percent)

# tso.exoplanet_spectrum.spectrum.wavelength_unit = u.micrometer
mask_use = tso.exoplanet_spectrum.spectrum.wavelength.mask

rebinned_wave = \
    tso.exoplanet_spectrum.spectrum.wavelength.data.value[~mask_use][4:-3:6]
rebinned_spec, rebinned_error = \
    spectres(rebinned_wave,
             tso.exoplanet_spectrum.spectrum.wavelength.data.value[~mask_use],
             tso.exoplanet_spectrum.spectrum.data.data.value[~mask_use],
             spec_errs=tso.exoplanet_spectrum.
             spectrum.uncertainty.data.value[~mask_use])

corrected_rebinned_spec, corrected_rebinned_error = \
    spectres(rebinned_wave,
             tso.exoplanet_spectrum.corrected_spectrum.wavelength.
             data.value[~mask_use],
             tso.exoplanet_spectrum.corrected_spectrum.data.
             data.value[~mask_use],
             spec_errs=tso.exoplanet_spectrum.corrected_spectrum.uncertainty.
             data.value[~mask_use])

from scipy import interpolate
f = interpolate.interp1d(wave_correction.value, corection_signal,
                         fill_value='extrapolate')
corrected_rebinned_spec = corrected_rebinned_spec*f(rebinned_wave)/1.03
corrected_rebinned_error = corrected_rebinned_error*f(rebinned_wave)/1.03
rebinned_spec = rebinned_spec*f(rebinned_wave)/1.03
rebinned_error = rebinned_error*f(rebinned_wave)/1.03

fig, ax = plt.subplots(figsize=(7, 5))
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
ax.plot(rebinned_wave, rebinned_spec, color="black", lw=5, alpha=0.9,
        zorder=5)
ax.errorbar(rebinned_wave, rebinned_spec, yerr=rebinned_error,
            fmt=".k", color="black", lw=5,
            alpha=0.9, ecolor="black",
            markeredgecolor='black', fillstyle='full', markersize=10,
            markerfacecolor='black', zorder=5)
ax.plot(wave_correction, corrected_signal_stevenson2014,
        color="red", lw=5, alpha=0.9, zorder=4)
ax.errorbar(wave_correction.value, corrected_signal_stevenson2014,
            yerr=error_corrected_signal_stevenson2014,
            fmt=".k", color="r", lw=5,
            alpha=0.9, ecolor="r",
            markeredgecolor='r', fillstyle='full', markersize=10,
            markerfacecolor='r', zorder=4)
ax.plot(wave_mandell, flux_mandell, color="r", lw=5, alpha=0.9, zorder=4)
#ax.errorbar(wave_mandell.value, flux_mandell.value, yerr=error_mandell.value,
#            fmt=".k", color="r", lw=5,
#            alpha=0.9, ecolor="r",
#            markeredgecolor='r', fillstyle='full', markersize=10,
#            markerfacecolor='r', zorder=4)
#ax.plot(wave_tsiaras, flux_tsiaras, color="y", lw=5, alpha=0.9, zorder=3)
#ax.errorbar(wave_tsiaras.value, flux_tsiaras.value, yerr=error_tsiaras.value,
#            fmt=".k", color="y", lw=5,
#            alpha=0.9, ecolor="y",
#            markeredgecolor='y', fillstyle='full', markersize=10,
#            markerfacecolor='y', zorder=3)
# ax.plot(tso.exoplanet_spectrum.spectrum.wavelength,
#        tso.exoplanet_spectrum.spectrum.data, lw=3, alpha=0.5, color='gray',
#        zorder=2)
# ax.errorbar(tso.exoplanet_spectrum.spectrum.wavelength.data.value,
#            tso.exoplanet_spectrum.spectrum.data.data.value,
#            yerr=tso.exoplanet_spectrum.spectrum.uncertainty.data.value,
#            fmt=".k", color='gray', lw=3, alpha=0.5, ecolor='gray',
#            markerfacecolor='gray',
#            markeredgecolor='gray', fillstyle='full', markersize=10,
#            zorder=2)
ax.plot(rebinned_wave, corrected_rebinned_spec, color="green", lw=4,
        alpha=0.9, zorder=6)
ax.errorbar(rebinned_wave, corrected_rebinned_spec,
            yerr=corrected_rebinned_error,
            fmt=".k", color="green", lw=4,
            alpha=0.9, ecolor="green",
            markeredgecolor='green', fillstyle='full', markersize=10,
            markerfacecolor='green', zorder=6)
axes = plt.gca()
axes.set_xlim([1.1, 1.7])
axes.set_ylim([1.25, 1.55])
ax.set_ylabel('Transit depth')
ax.set_xlabel('Wavelength')
plt.show()

residual_time_series = tso.calibration_results.residual.copy()
residual_unit = residual_time_series.data.unit
image0 = np.ma.array(residual_time_series.data,
                     mask=residual_time_series.mask)
image0.set_fill_value(np.nan)
image0 = image0[:, 0, :].filled().value
# model = tso.calibration_results.model_time_series[:,0,:].data.value
# image0 = image0/np.sqrt(model)
wave0 = tso.observation.dataset._wavelength
wave0_min = np.min(wave0)
wave0_max = np.max(wave0)
plt.imshow(image0,
           origin='lower',
           cmap='hot',
           interpolation='nearest',
           aspect='auto',
           extent=[0, image0.shape[1], wave0_min, wave0_max])
plt.colorbar().set_label("Residual ({})".format(residual_unit))
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
ax = autocorrelation_plot(df[10], lw=3, ax=ax, label="50")
ax = autocorrelation_plot(df[50], lw=3, ax=ax, label="150")
ax = autocorrelation_plot(df[90], lw=3, ax=ax, label="250")
ax = autocorrelation_plot(df, lw=3, ax=ax, label="All")
ax.legend(loc='best')
ax.set_title('ACF in wavelength direction')
plt.show()


nintegrations = image0.shape[0]
nwave = image0.shape[0]
fig, ax = plt.subplots(figsize=(7, 5))
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
s = pd.Series(np.random.normal(0.0, 1.0, nintegrations))
ax = autocorrelation_plot(s, lw=3, ax=ax, label="run1")
s = pd.Series(np.random.normal(0.0, 1.0, nintegrations))
ax = autocorrelation_plot(s, lw=3, ax=ax, label="run2")
s = pd.Series(np.random.normal(0.0, 1.0, nintegrations))
ax = autocorrelation_plot(s, lw=3, ax=ax, label="run3")
s = pd.Series(np.random.normal(0.0, 1.0, nintegrations))
ax = autocorrelation_plot(s, lw=3, ax=ax, label="run4")
s = pd.DataFrame(np.random.normal(0.0, 1.0, nintegrations*nwave).
                 reshape(nintegrations, nwave))
ax = autocorrelation_plot(s, lw=3, ax=ax, label="All")
ax.legend(loc='best')
ax.set_title('ACF in time diraction for random distributions')
plt.show()

fig, ax = plt.subplots(figsize=(7, 5))
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
s = pd.Series(np.random.normal(0.0, 1.0, nwave))
ax = autocorrelation_plot(s, lw=3, ax=ax, label="run1")
s = pd.Series(np.random.normal(0.0, 1.0, nwave))
ax = autocorrelation_plot(s, lw=3, ax=ax, label="run2")
s = pd.Series(np.random.normal(0.0, 1.0, nwave))
ax = autocorrelation_plot(s, lw=3, ax=ax, label="run3")
s = pd.Series(np.random.normal(0.0, 1.0, nwave))
ax = autocorrelation_plot(s, lw=3, ax=ax, label="run4")
s = pd.DataFrame(np.random.normal(0.0, 1.0, nintegrations*nwave).
                 reshape(nwave, nintegrations))
ax = autocorrelation_plot(s, lw=3, ax=ax, label="All")
ax.legend(loc='best')
ax.set_title('ACF in wavelength direction for random distributions')
plt.show()

# save planetary signal
tso.execute("save_results")

# plot planetary signal
tso.execute("plot_results")
