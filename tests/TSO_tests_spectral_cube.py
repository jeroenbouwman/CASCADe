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
                           "cascade_test_data_spectral_detector_cube.ini",
                           path=path)
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
            "cascade_test_data_spectral_detector_cube.ini", path=path)
print(tso.cascade_parameters)
print(cascade.initialize.cascade_configuration)
print(tso.cascade_parameters.isInitialized)

# load data
tso.execute("load_data")
plt.imshow(tso.observation.dataset.data[:, :, 15, 0])
plt.show()
plt.imshow(tso.observation.dataset_background.data[:, :, 15,  0])
plt.show()

plt.imshow(tso.observation.dataset._wavelength)
plt.show()
plt.imshow(tso.observation.dataset.wavelength[:, :, 15, 0])
plt.show()

with quantity_support():
    plt.plot(tso.observation.dataset.time.data[80, 18, 15, :],
             tso.observation.dataset.data.data[80, 18, 15, :])
    plt.show()

with quantity_support():
    plt.plot(tso.observation.dataset.wavelength.data[:, 18, 15, 0],
             tso.observation.dataset.data.data[:, 18, 15, 0])
    plt.show()

print(tso.observation.dataset.time.data.shape)
print(type(tso.observation.dataset.time.data))
print(tso.observation.dataset.time.data.unit)

# subtract background
tso.execute("subtract_background")
plt.imshow(tso.observation.dataset.data[:, :, 15, 0])
plt.show()

plt.imshow((tso.observation.dataset.mask.any(axis=2)).any(axis=2))
plt.show()

with quantity_support():
    plt.plot(tso.observation.dataset.time[80, 18, 15, :],
             tso.observation.dataset.data[80, 18, 15, :])
    plt.show()

# sigma clip data
tso.execute("sigma_clip_data")
plt.imshow(np.ma.median(tso.observation.dataset.data[:, :, 15, :], axis=2))
plt.show()
plt.imshow(tso.observation.dataset.mask[:, :, 15, :].all(axis=2))
plt.show()
with quantity_support():
    plt.plot(tso.observation.dataset.time[80, 18, 15, :],
             tso.observation.dataset.data[80, 18, 15, :])
    plt.show()

# eclipse model
tso.execute("define_eclipse_model")
plt.imshow(tso.model.light_curve_interpolated[:, 0, 15, :])
plt.show()
with quantity_support():
    plt.plot(tso.observation.dataset.time.data[80, 18, 15, :],
             tso.model.light_curve_interpolated[80, 18, 15, :])
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

plt.plot(tso.cpm.spectral_trace)
plt.show()
plt.plot(tso.cpm.position)
plt.show()

plt.plot(tso.observation.spectral_trace['wavelength_pixel'].value,
         tso.observation.spectral_trace['positional_pixel'].value -
         tso.cpm.spectral_trace)
axes = plt.gca()
axes.set_ylim([-1, 1])
plt.show()

