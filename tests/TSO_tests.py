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
import warnings
from sklearn.preprocessing import RobustScaler
from scipy.linalg import pinv2

warnings.simplefilter("ignore")

sns.set_style("white")

# create transit spectoscopy object
tso = cascade.TSO.TSOSuite()

# initialization with ini files
path = cascade.initialize.default_initialization_path
tso = cascade.TSO.TSOSuite("cascade_test_cpm_spectra.ini",
                           "cascade_test_object.ini",
                           "cascade_test_data_spectra.ini", path=path)

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

# create TSO object and initialize
tso = cascade.TSO.TSOSuite()
path = cascade.initialize.default_initialization_path
tso.execute("initialize", "cascade_test_cpm_spectra.ini",
            "cascade_test_object.ini",
            "cascade_test_data_spectra.ini", path=path)

print(tso.cascade_parameters)
print(cascade.initialize.cascade_configuration)
print(tso.cascade_parameters.isInitialized)
assert tso.cascade_parameters == cascade.initialize.cascade_configuration
assert tso.cascade_parameters.isInitialized is True

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

with quantity_support():
    plt.plot(tso.observation.dataset.time.data[80, :],
             tso.observation.dataset.data.data[80, :])
    plt.show()

with quantity_support():
    plt.plot(tso.observation.dataset.wavelength.data[:, 0],
             tso.observation.dataset.data.data[:, 0])
    plt.show()

print(tso.observation.dataset.time.data.shape)
print(type(tso.observation.dataset.time.data))
print(tso.observation.dataset.time.data.unit)
assert tso.observation.dataset.time.data.shape == (130, 280)
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

plt.imshow(image_test_roi,
           origin='lower',
           cmap='hot',
           interpolation='nearest',
           aspect='auto')
plt.show()
plt.imshow(cleaned_image_test_roi.data.value,
           origin='lower',
           cmap='hot',
           interpolation='nearest',
           aspect='auto')
plt.show()

plt.imshow(RobustScaler().fit_transform(cleaned_image_test_roi),
           origin='lower',
           cmap='hot',
           interpolation='nearest',
           aspect='auto')
plt.colorbar().set_label("Intensity")
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

from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans
from astropy.convolution import Gaussian1DKernel
kernel = tso.observation.instrument_calibration.convolution_kernel
kernel_size = kernel.shape[0]
kernel_1d = Gaussian1DKernel(0.2, x_size=kernel_size)
kernel = np.repeat(np.expand_dims(kernel, axis=1),
                    (kernel_size), axis=1)
kernel = kernel*kernel_1d.array[:, None].T
kernel = kernel/np.sum(kernel)

image0b = np.ma.array(tso.observation.dataset.data.data.value.copy(),
                      mask=tso.observation.dataset.mask)
image0b.set_fill_value(np.nan)
RS = RobustScaler(with_scaling=True)
RS.fit(image0b.T)
X_scaled_masked = RS.transform(image0b.T)
bla = np.ma.array(X_scaled_masked, mask=image0b.T.mask)
bla.set_fill_value(np.nan)
bla2 = interpolate_replace_nans(bla.filled().T, kernel, boundary='extend')
plt.plot(bla2.T)
plt.plot(np.ma.median(bla2.T, axis=1))
plt.ylim([-4, 4])
plt.show()

################################################################

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

print(tso.model.transit_timing)
assert tso.model.transit_timing == [0.4827232723272327, 0.5172767276727672]

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
    mask_use = ~tso.observation.dataset.data[80, :].mask
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(tso.observation.dataset.time.data[80, mask_use],
            tso.observation.dataset.data.data[80, mask_use], label='Data')
    ax.scatter(tso.observation.dataset.time.data[80, mask_use],
            tso.calibration_results.model_time_series.data[80, 0, mask_use],
            label='Model')
    ax.legend(loc='best')
    ax.set_title('Comparison Data and Regression Model')
    ax.set_xlabel('Phase')
    plt.show()

# extract planetary signal
tso.execute("extract_spectrum")

fig, ax = plt.subplots(figsize=(7, 5))
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
ax.plot(tso.exoplanet_spectrum.spectrum.wavelength,
        np.ma.abs(tso.exoplanet_spectrum.weighted_image), lw=3)
ax.set_ylabel('| Signal * weight |')
ax.set_xlabel('Wavelength')
plt.show()

# correct the extracted planetary signal for non uniform subtraction
# of averige signal
tso.execute("correct_extracted_spectrum")

################################################
# Read resuls from literature/other methods
################################################
wavelength_corrected_spectrum = (tso.exoplanet_spectrum.corrected_spectrum.
                                 wavelength)
corrected_spectrum = (tso.exoplanet_spectrum.corrected_spectrum.data)

error_corrected_spectrum = (tso.exoplanet_spectrum.corrected_spectrum.
                            uncertainty)

path_old = '/home/bouwman/SST_OBSERVATIONS/projects_HD189733/REDUCED_DATA/'
spec_instr_model_ian = ascii.read(path_old+'results_ian.dat', data_start=1)
wave_ian = (spec_instr_model_ian['lam_micron']*u.micron)
flux_ian = (spec_instr_model_ian['depthA'] *
            u.dimensionless_unscaled).to(u.percent)
error_ian = (spec_instr_model_ian['errorA'] *
             u.dimensionless_unscaled).to(u.percent)

fig, ax = plt.subplots(figsize=(10, 7))
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
ax.plot(wavelength_corrected_spectrum,
        corrected_spectrum,
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

residual_time_series = tso.calibration_results.residual.copy()
residual_unit = residual_time_series.data.unit
image0 = np.ma.array(residual_time_series.data,
                     mask=residual_time_series.mask)
image0.set_fill_value(np.nan)
image0 = image0[:, 0, :].filled().value
wave0 = np.ma.median(tso.observation.dataset.wavelength, axis=1)
iw_min = np.ma.argmin(wave0)
iw_max = np.ma.argmax(wave0)
wave0_min = wave0[iw_min].value
wave0_max = wave0[iw_max].value
plt.imshow(image0[iw_min:iw_max, :],
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
s = pd.DataFrame(np.random.normal(0.0, 1.0, 280*128).reshape(280, 128))
ax = autocorrelation_plot(s, lw=3, ax=ax, label="All")
ax.legend(loc='best')
ax.set_title('ACF in time diraction for random distributions')
plt.show()

fig, ax = plt.subplots(figsize=(7, 5))
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
s = pd.Series(np.random.normal(0.0, 1.0, 128))
ax = autocorrelation_plot(s, lw=3, ax=ax, label="run1")
s = pd.Series(np.random.normal(0.0, 1.0, 128))
ax = autocorrelation_plot(s, lw=3, ax=ax, label="run2")
s = pd.Series(np.random.normal(0.0, 1.0, 128))
ax = autocorrelation_plot(s, lw=3, ax=ax, label="run3")
s = pd.Series(np.random.normal(0.0, 1.0, 128))
ax = autocorrelation_plot(s, lw=3, ax=ax, label="run4")
s = pd.DataFrame(np.random.normal(0.0, 1.0, 280*128).reshape(128, 280))
ax = autocorrelation_plot(s, lw=3, ax=ax, label="All")
ax.legend(loc='best')
ax.set_title('ACF in wavelength direction for random distributions')
plt.show()

# save planetary signal
tso.execute("save_results")

# plot planetary signal
tso.execute("plot_results")

#########################################################

cal_file = 'IRSX_SL1_S18.18.0_cal_rsrf_default_droop.tbl'
path_cal = '/home/bouwman/IRS_Calibration/SPITZER_IA_CAL/S18.18.0/'
rsrf = ascii.read(path_cal+cal_file, data_start=0)
wave_rsrf = rsrf['col1'].data
rsrf_1 = rsrf['col2'].data
rsrf_2 = rsrf['col4'].data

wave = np.ma.array(tso.exoplanet_spectrum.spectrum.wavelength.data,
                   mask=tso.exoplanet_spectrum.spectrum.wavelength.mask)
star = 0.5 * (wave/7.45)**-2
rsrf_new = (tso.calibration_results.parameters[-1, :, 0]) / \
    (tso.calibration_results.signal[:, 0])
rsrf_old = (rsrf_1 + rsrf_2)*0.47*star
plt.plot(wave, rsrf_new)
plt.plot(wave, rsrf_old)
plt.show()

plt.plot(wave, rsrf_new/rsrf_old)
plt.show()

#####################################################
