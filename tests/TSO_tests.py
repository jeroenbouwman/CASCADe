# -*- coding: utf-8 -*-

import cascade
from matplotlib import pyplot as plt
from astropy.visualization import quantity_support
import numpy as np

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

# reset parameters
tso.execute("reset")
print(tso.cascade_parameters.isInitialized)

# initialize with providing ini files
tso.execute("initialize")
print(tso.cascade_parameters.isInitialized)

# create TSO object and initialize
tso = cascade.TSO.TSOSuite()
path = cascade.initialize.default_initialization_path
tso.execute("initialize", "cascade_test_cpm.ini", "cascade_test_object.ini",
            "cascade_test_data_spectra.ini", path=path)
print(tso.cascade_parameters)
print(cascade.initialize.cascade_configuration)
print(tso.cascade_parameters.isInitialized)

# load data
tso.execute("load_data")
plt.imshow(tso.observation.dataset.data[:, :])
plt.show()


plt.imshow(tso.observation.dataset._wavelength)
plt.show()
plt.imshow(tso.observation.dataset.wavelength[:, :])
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

# subtract background
tso.execute("subtract_background")
plt.imshow(tso.observation.dataset.data[:, :])
plt.show()

plt.imshow(tso.observation.dataset.mask)
plt.show()

with quantity_support():
    plt.plot(tso.observation.dataset.time[80, :],
             tso.observation.dataset.data[80, :])
    plt.show()

# sigma clip data
tso.execute("sigma_clip_data")
plt.imshow(tso.observation.dataset.data[:, :])
plt.show()
plt.imshow(tso.observation.dataset.mask[:, :])
plt.show()
with quantity_support():
    plt.plot(tso.observation.dataset.time[80, :],
             tso.observation.dataset.data[80, :])
    plt.show()

# eclipse model
tso.execute("define_eclipse_model")
plt.imshow(tso.model.light_curve_interpolated[:, :])
plt.show()
with quantity_support():
    plt.plot(tso.observation.dataset.time.data[80, :],
             tso.model.light_curve_interpolated[80, :])
    plt.show()

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
print(len(tso.cpm.regressor_list))
print(len(tso.cpm.regressor_list[0]))
print(len(tso.cpm.regressor_list[0][0]))
print(tso.cpm.regressor_list[0][0])
print(tso.cpm.regressor_list[0][0][0])
print(tso.cpm.regressor_list[0][0][1])


# setup of regression matrix
tso.return_all_design_matrices()
print(tso.cpm.design_matrix[0][0][0].shape)
mask = tso.cpm.design_matrix[0][50][0].mask
plt.imshow(mask)
plt.show()
data_masked = tso.cpm.design_matrix[0][50][0]
plt.imshow(data_masked)
plt.show()

tso.return_all_design_matrices(clip=True)
print(tso.cpm.design_matrix[0][0][0].shape)
mask = tso.cpm.design_matrix[0][50][0].mask
plt.imshow(mask)
plt.show()
data_masked = tso.cpm.design_matrix[0][50][0]
plt.imshow(data_masked)
plt.show()


# create calibrated time series and derive planetary signal
tso.execute("calibrate_timeseries")
print(np.ma.max(tso.calibration_results.regularization))
print(np.ma.min(tso.calibration_results.regularization))
plt.imshow(tso.calibration_results.regularization)
plt.show()

extent = [0, 127, 15.5, 7.5]
delta_lam = 15.5-7.5
image = (tso.calibration_results.signal/
         tso.calibration_results.error_signal**2).copy()
fig, ax = plt.subplots(figsize=(6, 6))
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
#plt.cm.Reds
cmap = plt.cm.gist_heat
cmap.set_bad('black', 1.)
ax.imshow(np.ma.abs(image), origin='lower', aspect='auto',
          cmap=cmap, interpolation='none', extent=extent, vmin=0, vmax=1000)
ax.set_aspect(128.0 / delta_lam)
ax.set_xlabel('Pixel Number Spatial Direction')
ax.set_ylabel('Wavelength')

###################################################
# calculate median spectrum over imaging array    ##
# to extract the spectral signature of the planet ##
###################################################
wave_use_image = tso.observation.dataset.wavelength.data.value[:,0]
extraction_mask = tso.cpm.extraction_mask[0]
masked_wave_use_image = np.ma.array(wave_use_image, mask=extraction_mask)
#masked_wave_use_image = np.ma.array(wave_use_image.data.value, mask=wave_use_image.mask)
masked_wave_use_image = masked_wave_use_image[:,None]

npix,mpix = masked_wave_use_image.shape
av = np.ma.average(tso.calibration_results.signal, axis=1,
                   weights=(np.ma.ones((npix, mpix)) /
                            tso.calibration_results.error_signal)**2)
av_error = np.ma.ones((npix)) / np.ma.sum((np.ma.ones((npix, mpix)) /
                                          tso.calibration_results.error_signal)**2, axis=1)
av_error = np.ma.sqrt(av_error)
av_wave = np.ma.average(masked_wave_use_image, axis=1,
                        weights=(np.ma.ones((npix, mpix)) /
                                 tso.calibration_results.error_signal)**2)

median_eclipse_depth = 0.00435   # SL1
spectrum = (av * (1.0 + median_eclipse_depth) +
            median_eclipse_depth)
error_spectrum = av_error * (1.0 + median_eclipse_depth)

path_old = '/home/bouwman/SST_OBSERVATIONS/projects_HD189733/REDUCED_DATA/'
spec_instr_model_ian = ascii.read(path_old+'results_ian.dat', data_start=1)
fig, ax = plt.subplots(figsize=(10, 6))
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
ax.plot(spec_instr_model_ian['lam_micron'], spec_instr_model_ian['depthA'],
            color="r", lw=2, alpha=0.8)
ax.errorbar(spec_instr_model_ian['lam_micron'],
            spec_instr_model_ian['depthA'],
            yerr=spec_instr_model_ian['errorA'],
            fmt=".k", color="r", lw=2,
            alpha=0.8, ecolor="r")

ax.plot(av_wave, spectrum, lw=3, alpha=0.7, color='blue')
ax.errorbar(av_wave, spectrum, yerr=error_spectrum, fmt=".k", color='blue',
            lw=3, alpha=0.7, ecolor='blue', mfc='blue')

axes = plt.gca()
axes.set_xlim([7.5, 15.5])
axes.set_ylim([-0.00, 0.008])
ax.set_ylabel('Fp/Fstar')
ax.set_xlabel('Wavelength')

