# -*- coding: utf-8 -*-
import cascade
from matplotlib import pyplot as plt

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
# test on general interface level
observation = cascade.instruments.Observation()
plt.imshow(observation.data.data[:, :, 0])
plt.show()
