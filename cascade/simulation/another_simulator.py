# another simulator


import numpy as np
import os
import matplotlib.pyplot as plt
import astropy.units as u
from configobj import ConfigObj
import batman
import ast
from astropy.io import fits
from numpy.random import poisson
import astropy.units as u


####
from cascade.utilities.readwrite_spectra import plot_spectrum_2d
from cascade.utilities.readwrite_spectra import read_depth_ecsv
from cascade.utilities.readwrite_spectra import write_depth_ecsv
from cascade.utilities.readwrite_spectra import write_spectra_dir

from cascade.utilities.noises import add_noise
#from cascade.utilities.noises import add_sine_noise
from read_fluxes import read_fluxes
from compute_star_flux_evolution import compute_star_flux_evolution
from compute_planet_flux_evolution import compute_planet_flux_evolution


def return_par_from_ini_ter(file_object, file_model, verbose=False):
    """
    Get relevant parameters for flux of the planet from CASCADe
    intitialization files
    
    Parameters
    ----------
    file_object : 'str'
    name of the object configuration file,
    file_model : 'str'
    name of the model configuration file
    
    Returns
    -------
    :param quantity planet_radius : radius of the planet, unit : a length unit, e.g. jupiter radius
    :param quantity planet_semi_axis : semi-major axis of the planet, unit : a length unit, e.g. astronomical unit
    :param float albedo : albedo, number between zero (no reflection) and one
    """

    configM = ConfigObj(file_model)
    pattern = configM['OBSERVATIONS']['observations_id']
    save_path = configM['OBSERVATIONS']['observations_path']
    target_name = configM['OBSERVATIONS']['observations_target_name']
    fluxes_file = configM['SIMULATION']['fluxes_file']
    observations_id = configM['OBSERVATIONS']['observations_id']

    #######
    configO = ConfigObj(file_object)
    if(verbose): print('config object', configO['OBJECT'])
    planet_period = u.Quantity(configO['OBJECT']['object_period'])
    planet_ephemeris = u.Quantity(configO['OBJECT']['object_ephemeris'])
    return planet_period, planet_ephemeris, save_path, target_name, fluxes_file, observations_id


### compute of the phase
def compute_phase(file_object, file_model, verbose=False):
    configM = ConfigObj(file_model)
    model_phase_range = float(configM['MODEL']['model_phase_range'])
    model_nphase_points = int(configM['MODEL']['model_nphase_points'])
    sampling_time = u.Quantity(configM['SIMULATION']['sampling_time'])
    #
    configO = ConfigObj(file_object)
    planet_period = u.Quantity(configO['OBJECT']['object_period'])
    planet_ephemeris = u.Quantity(configO['OBJECT']['object_ephemeris'])
    #
    model_ntime_point = (model_phase_range*planet_period/sampling_time).decompose().value
    if(model_ntime_point > model_nphase_points/2):print('problem with model_nphase_points')
    model_ntime_point = int(model_ntime_point)
    phase = np.linspace(-model_phase_range/2, model_phase_range/2, model_ntime_point)
    #
    return phase, sampling_time

def another_simulator(file_object, file_model, verbose=False, plot=False):

    planet_period, planet_ephemeris, save_path, target_name, fluxes_file, observations_id = return_par_from_ini_ter(file_object, file_model, verbose=False)
    
    # to be replace by the actual comput of the fluxes with BlackBody approximation
    #######  Read the fluxes of the planet and the star  1D (wavelenght)  #######
    wavelength, star_flux1d, planet_day_flux1d, planet_night_flux1d = read_fluxes(fluxes_file)
    print(wavelength.shape, star_flux1d.shape, planet_day_flux1d.shape, planet_night_flux1d.shape)

    #  compute the phase, the time dependancy of the fluxes
    phase, sampling_time = compute_phase(file_object, file_model, verbose=False)

    ######  compute the flux of the planet  2D (wavelenght+ time) ####################
    planet_flux2d = compute_planet_flux_evolution(star_flux1d,  planet_day_flux1d, planet_night_flux1d,
                      phase, file_object, file_model, verbose=verbose)
    if(plot):plot_spectrum_2d(planet_flux2d.value, wavelength, phase, 'planet flux', normalise=False)

    #########  compute the flux of the star 2D (wavelenght+time)############
    star_flux2d = compute_star_flux_evolution(star_flux1d, wavelength, phase, file_object, file_model)
    if(plot):plot_spectrum_2d(star_flux2d.value, wavelength, phase, 'star flux', normalise=False)

    #######  Total Flux ###########
    total_flux2d = star_flux2d + planet_flux2d
    plot_spectrum_2d(total_flux2d.value, wavelength, phase, 'total flux', normalise=False)

    ### add Poisson noise ###########
    total_flux2de = total_flux2d*sampling_time
    noisy_flux2de = add_noise(total_flux2de.value)*total_flux2de.unit
    plot_spectrum_2d(noisy_flux2de.value, wavelength, phase, 'noisy flux electons', normalise=False)


    #########  Write the output #######
    # time is a reserved function
    times = phase*planet_period+planet_ephemeris
    out_dir = os.path.join(save_path, target_name)
    out_dir = os.path.join(out_dir,'SPECTRA')
    os.makedirs(out_dir,  exist_ok=True)
    uncertainties = np.sqrt(total_flux2de.value)*total_flux2de.unit
    mask = np.full((total_flux2de.shape), False, dtype=bool)
    pattern=observations_id
    write_spectra_dir(out_dir, pattern, wavelength, noisy_flux2de, uncertainties, mask, times, verbose=False, overwrite=True)
    return

