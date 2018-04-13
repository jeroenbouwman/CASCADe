# -*- coding: utf-8 -*-

import cascade
from matplotlib import pyplot as plt
from astropy.visualization import quantity_support
import numpy as np
from astropy.io import ascii
from astropy import units as u
from astropy.units import cds
from scipy import stats


# create transit spectoscopy object
tso = cascade.TSO.TSOSuite()

# initialization with ini files
path = cascade.initialize.default_initialization_path
tso = cascade.TSO.TSOSuite("cascade_test_cpm.ini",
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
tso.execute("initialize", "cascade_test_cpm.ini", "cascade_test_object.ini",
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
           extent =[0, image0.shape[1], wave0_min, wave0_max])
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
           extent =[0, image0.shape[1], wave0_min, wave0_max])
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
           extent =[0, image0.shape[1], wave0_min, wave0_max])
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

# eclipse model
tso.execute("define_eclipse_model")

plt.imshow(tso.model.light_curve_interpolated[0][:, :],
           origin='lower',
           cmap='hot',
           interpolation='nearest',
           aspect='auto',
           extent =[0, image0.shape[1], wave0_min, wave0_max])
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

with quantity_support():
    plt.plot(tso.observation.spectral_trace['wavelength_pixel'],
             tso.observation.spectral_trace['positional_pixel'])
    plt.show()
with quantity_support():
    plt.plot(tso.observation.spectral_trace['wavelength_pixel'],
             tso.observation.spectral_trace['wavelength'])
    plt.show()
with quantity_support():
    plt.plot(tso.observation.dataset.time[80, :],
             tso.observation.dataset.position[80, :])
    plt.show()
plt.plot(tso.cpm.spectral_trace)
plt.show()
plt.plot(tso.cpm.position[80, :])
plt.show()

# set the extraction area
tso.execute("set_extraction_mask")

print(tso.cpm.extraction_mask[0].shape)
plt.plot(tso.cpm.extraction_mask[0])
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
#       second: [# data pixel in wavelength direction]
#       third: [0=pixel coord; 1=list of regressors]
#       forth: [0=coordinate wave direction; 1=coordinate spatial direction]
#   design_matrix:
#       first: [# nod]
#       second: [# data pixel in wavelength direction]
#       third: [0]

i_data_point_in_wavelength_direction = 50
reg_list = \
    tso.cpm.regressor_list[0][i_data_point_in_wavelength_direction][1][0]
reg_matrix = tso.cpm.design_matrix[0][i_data_point_in_wavelength_direction][0]

plt.imshow(reg_matrix.mask)
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
tso.return_all_design_matrices(clip=True, clip_pctl_time=0.08,
                               clip_pctl_regressors=0.04)

reg_matrix = tso.cpm.design_matrix[0][i_data_point_in_wavelength_direction][0]

plt.imshow(reg_matrix.mask)
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

#for i in range(71):
#    ks, prob = stats.ks_2samp(flux_ian.value,
#                          (tso.exoplanet_spectrum.spectrum.data.data.
#                           value[~tso.exoplanet_spectrum.spectrum.mask])[::-1] -
#                           (0.035-i*0.001))
#    print(-(0.035-i*0.001), ks, prob)


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

wave = tso.exoplanet_spectrum.spectrum.wavelength.data.value
star = 0.5 * (wave/wave[0])**-2
rsrf_new = (tso.calibration_results.parameters[0, :, 0])/(tso.calibration_results.signal[:,0])
rsrf_old = (rsrf_1 + rsrf_2)*0.47*star
plt.plot(wave, rsrf_new)
plt.plot(wave, rsrf_old)
plt.show()

plt.plot(wave, rsrf_new/rsrf_old)
plt.show()



#####################################################

cal_signal_depth = float(tso.cascade_parameters.cpm_calibration_signal_depth)
mean_eclipse_depth = float(tso.cascade_parameters.observations_median_signal)

plt.plot(tso.calibration_results.calibrated_time_series[50, :].T)
plt.show()

plt.plot(tso.exoplanet_spectrum.wavelength,
         tso.calibration_results.calibration_signal *
         (1+cal_signal_depth)+cal_signal_depth)
axes = plt.gca()
axes.set_ylim([cal_signal_depth*0.0, cal_signal_depth*5])
plt.show()

plt.plot(tso.exoplanet_spectrum.wavelength,
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
