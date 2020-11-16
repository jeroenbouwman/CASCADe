#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# This file is part of the CASCADe package which has been
# developed within the ExoplANETS-A H2020 program.

##  It uses as input the spectra files produced by simulate_spectra_correl.py
##  It can check the impact of the parameter cpm_relative_sig_value_limit on the reconstruction
#
import time
import matplotlib.pyplot as plt
import os
import numpy as np
from astropy.io import ascii
import astropy.units as u
#
file_model='init_files/cascade_WASP43b_calibrate_Atlasstar.ini'
file_object='init_files/cascade_WASP43b_object.ini'
#
os.environ["CASCADE_WARNINGS"] = 'on'
os.environ["CASCADE_INITIALIZATION_FILE_PATH"] = os.getcwd()
os.environ['CASCADE_DATA_PATH'] = os.getcwd()
SAVE = True
#
import cascade
#
def ameliore( tso, rcond_limit, reference):
    tso.cascade_parameters.cpm_relative_sig_value_limit = rcond_limit
    tso.execute("correct_extracted_spectrum")
    diff = tso.exoplanet_spectrum.corrected_spectrum.data.data.value - reference
    #print("mean={:.2e}, standard deviation={:.2e}, rcond_limit={:.2e}".format(diff.mean(), diff.std(), rcond_limit))
    return diff.std()

#################  BEGINS HERE

#init_path = os.path.join(os.path.dirname(cascade.__path__[0]), 'examples/init_files')
#init_path = os.getcwd()

# create transit spectroscopy object
tso = cascade.TSO.TSOSuite()
tso.execute("reset")

tso.execute("initialize", file_object, file_model)#, path=init_path)

print('observation path', tso.cascade_parameters.observations_path)
print('save path', tso.cascade_parameters.cascade_save_path)

## read the reference transit depth
atmosphere_file = tso.cascade_parameters.atmosphere_file
print('atmosphere_file', atmosphere_file)
data = ascii.read(atmosphere_file)
plt.figure()
plt.plot(data['wavelength'], data['depth'])
plt.xlabel("wavelength")
plt.ylabel("transit depth")

# load the spectral data
tso.execute("load_data")

tso.execute("subtract_background")

# filter data and create cleaned dataset
tso.execute("filter_dataset")

# determine the relative movement from the input data.
# In case of 1D spectra, positional data given in the fits header will be used.
tso.execute("determine_source_movement")

# set the extraction area
tso.execute("set_extraction_mask")

# setup regressors
tso.execute("select_regressors")

# eclipse model
tso.execute("define_eclipse_model")

# create calibrated time series and derive planetary signal
tso.execute("calibrate_timeseries")

# extract planetary signal
tso.execute("extract_spectrum")

# correct the extracted planetary signal for non uniform
#  subtraction of average signal
tso.execute("correct_extracted_spectrum")
###################

target_name = tso.observation.dataset_parameters['obs_target_name']
# WASP43b
#rcond_limit= float(tso.cascade_parameters.cpm_relative_sig_value_limit)
#print("rcond_limit ={:.2e}".format( rcond_limit))
rcond_limit= tso.cascade_parameters.cpm_relative_sig_value_limit

####  compute the reference transit depth by the ratio of planet radius on star radius (squared)
planet_radius = u.Quantity(tso.cascade_parameters.object_radius)
#'1.036 Rjup'
star_radius = u.Quantity(tso.cascade_parameters.object_radius_host_star)
# '0.667 Rsun'
rp = (planet_radius/star_radius).decompose().value
ref_transit = rp*rp*100
# 2.547
reference = ref_transit*data['depth']**2

#### plot the transit depth, reconstructed transit depth and theoric transit depth
plt.style.use('classic')
plt.figure()
plt.title('result of cascade '+target_name+' rcond_limit='+rcond_limit)
plt.plot(tso.exoplanet_spectrum.corrected_spectrum.wavelength, tso.exoplanet_spectrum.corrected_spectrum.data,'-p',label='corrected spectrum')
plt.plot(tso.exoplanet_spectrum.spectrum.wavelength, tso.exoplanet_spectrum.spectrum.data,'o', label='spectrum')
plt.plot(data['wavelength'], reference, label='reference', linewidth=2)
plt.legend()
plt.xlabel("wavelength")
plt.ylabel("relative depth transit")
plt.show()
if (SAVE): plt.savefig('result_cascade.png')

ameliore( tso, rcond_limit, reference)
# mean=-1.95e-03, standard deviation=3.78e-02, rcond_limit=1.00e-03
num = 1000
x = np.linspace(1e-5, 1e-2, num=num)
num = 100
x = np.linspace(1e-2, 1, num=num)

print(x.min(), x.max())
y = np.zeros(num)
for i in np.arange(num):
    y[i] =  ameliore( tso, x[i], reference)
plt.figure()
plt.semilogx(x, y)
plt.semilogx(x, y, '+')
plt.title("correct_extracted_spectrum wasp49b")
plt.xlabel('cpm_relative_sig_value_limit')
plt.ylabel('standard deviation of the corrected spectrum')
plt.savefig('plot_rcond_limit_fake_wasp49b.png')


# save planetary signal
tso.execute("save_results")

# plot planetary signal
tso.execute("plot_results")


