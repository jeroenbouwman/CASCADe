# IPython log file

from compute_planet_flux_evolution import *

file_model='init_files/cascade_WASP43b_calibrate_Atlasstar.ini'
file_object='init_files/cascade_WASP43b_object.ini'
nw = 100
nt = 600
wavelength = (np.arange(nw)/10.+4.)*u.micron
star_flux = np.ones(nw)*u.electron/u.second
phase = np.arange(nt)/1000.-0.3
star_flux = np.ones(nw)*u.electron/u.second*100.
star_flux = np.ones(nw)*u.electron/u.second*1000.
planet_flux_day_0 = np.ones(nw)*u.electron/u.second*10.
planet_flux_night_0 = np.ones(nw)*u.electron/u.second*1.
planet_flux = compute_planet_flux_evolution(star_flux,  planet_flux_day_0, planet_flux_night_0,\
                          phase, file_object, file_model, verbose=False)
planet_flux.shap
planet_flux.shape
from cascade.utilities.readwrite_spectra import plot_spectrum_2d
plot_spectrum_2d(planet_flux.value, wavelength, phase, 'planet flux', normalise=False)
get_ipython().run_line_magic('logstart', 'check_compute_flux_star_evolution.py')
get_ipython().run_line_magic('logstop', '')
