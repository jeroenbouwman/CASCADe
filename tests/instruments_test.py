# -*- coding: utf-8 -*-
import cascade
from matplotlib import pyplot as plt
from astropy.visualization import quantity_support

# initialize cascade
path = cascade.initialize.default_initialization_path
cascade_param = \
    cascade.initialize.configurator(path+"cascade_test_cpm.ini",
                                    path+"cascade_test_object.ini",
                                    path+"cascade_test_data_spectral_images.ini")

# Test on observatory level
observatory = cascade.instruments.Spitzer()
plt.imshow(observatory.data.data[:, :, 0])
plt.show()
plt.imshow(observatory.data_background.data[:, :, 0])
plt.show()
with quantity_support():
    plt.plot(observatory.spectral_trace['wavelength_pixel'],
             observatory.spectral_trace['positional_pixel'])
    plt.show()

# test on general interface level
observation = cascade.instruments.Observation()
plt.imshow(observation.dataset.data[:, :, 0])
plt.show()
plt.imshow(observation.dataset_background.data[:, :, 0])
plt.show()
with quantity_support():
    plt.plot(observation.spectral_trace['wavelength_pixel'],
             observation.spectral_trace['positional_pixel'])
    plt.show()

# Tests for detector cube data
cascade_param.reset()
# initialize cascade
path = cascade.initialize.default_initialization_path
cascade_param = \
    cascade.initialize.configurator(path+"cascade_test_cpm.ini",
                                    path+"cascade_test_object.ini",
                                    path+"cascade_test_data_spectral_detector_cube.ini")
# test on general interface level
observation = cascade.instruments.Observation()
observation.dataset.data.shape

plt.imshow(observation.dataset.data[:, :, 15, 0])
plt.show()
plt.imshow(observation.dataset_background.data[:, :,15, 0])
plt.show()

with quantity_support():
    plt.plot(observation.dataset.time.data[80, 18, :, 0],
             observation.dataset.data.data[80, 18, :, 0])
    plt.show()

with quantity_support():
    plt.plot(observation.dataset.time.data[80, 18, 15, :],
             observation.dataset.data.data[80, 18, 15, :])
    plt.show()

