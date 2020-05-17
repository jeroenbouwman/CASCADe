# compute_star_flux_evolution
from cascade.utilities.readwrite_spectra import plot_spectrum_2d
from compute_star_flux_evolution import compute_star_flux_evolution
import numpy as np
import os
import matplotlib.pyplot as plt
import astropy.units as u
#

file_model='init_files/cascade_WASP43b_calibrate_Atlasstar.ini'
file_object='init_files/cascade_WASP43b_object.ini'
nw = 100
nt = 600
wavelength = (np.arange(nw)/10.+4.)*u.micron
star_flux = np.ones(nw)*u.electron/u.second
phase = np.arange(nt)/1000.-0.3

star_flux2d = compute_star_flux_evolution(star_flux, wavelength, phase, file_object, file_model)

plot_spectrum_2d(star_flux2d.value, wavelength, phase, 'star flux', normalise=False)
