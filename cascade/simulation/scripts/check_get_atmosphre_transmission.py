# check get atmosphe transmission

from compute_star_flux_evolution import *
import numpy as np
import os
import matplotlib.pyplot as plt
import astropy.units as u
#

file_model='init_files/cascade_WASP43b_calibrate_Atlasstar.ini'
file_object='init_files/cascade_WASP43b_object.ini'
wavelength = (np.arange(100)/10.+4.)*u.micron
#
depths = get_atmosphere_transmission(wavelength, file_model, plot=True, verbose=True )
