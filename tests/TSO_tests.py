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
tso.execute("setup_design_matrix")
print(tso.cpm.design_matrix[0][0][0].shape)
print(tso.cpm.design_matrix[0][1][1].shape)