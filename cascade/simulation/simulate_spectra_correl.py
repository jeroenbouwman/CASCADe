#  Create dummy light curves with systematics

#####  see simulate_spectra_correl
### Replace CML.transit_depth_create_fake_lc, create_fake_lc_causal_noise
### Replace  exoplanet_tools class batman_model  which neeeds the singleton class config_param

import numpy as np
import os
import matplotlib.pyplot as plt
import astropy.units as u
from configobj import ConfigObj
import batman
import ast
from astropy.io import fits
from numpy.random import poisson

####
from cascade.utilities.readwrite_spectra import plot_spectrum_2d
from cascade.utilities.readwrite_spectra import read_depth_ecsv
from cascade.utilities.readwrite_spectra import write_spectra_dir
from cascade.utilities.noises import add_noise
from cascade.utilities.noises import add_sine_noise

#######################

def return_par_from_ini(file_object='cascade_WASP43b_object.ini',
                        file_model='cascade_WASP43b_calibrate_planet_spectrum.ini'):
    """
    Get relevant parameters for lightcurve model from CASCADe
    intitialization files
    
    Parameters
    ----------
    file_object : 'str'
        optional, name of the object configuration file,
    file_model : 'str'
        optional, name of the model configuration file

    Returns
    -------
    params: 'batman parameter object'
        model parameters for batman lightcurve mode (object)
    phase: 'numpy array'
        Orbital phase of planet where the lightcurve should be computed
    planet_period: 'quantity' 
        period of the planet
    planet_ephemeris: 'quantity'
        date of the middle of one planet transit (Modified  Barycentered Julian Date)
    
    """
    config = ConfigObj(file_model)
    print('MODEL', config['MODEL'])
    model_phase_range = float(config['MODEL']['model_phase_range'])
    model_nphase_points = int(config['MODEL']['model_nphase_points'])
    model_limb_darkening=config['MODEL']['model_limb_darkening']
    model_limb_darkening_coeff=config['MODEL']['model_limb_darkening_coeff']
    chaine = model_limb_darkening_coeff[0]+','+model_limb_darkening_coeff[1]
    model_limb_darkening_coeff = ast.literal_eval(chaine)
    print(model_limb_darkening)
    #######
    config = ConfigObj(file_object)
    print('config object', config['OBJECT'])
    planet_name = config['OBJECT']['object_name']
    planet_radius = u.Quantity(config['OBJECT']['object_radius'])
    star_radius =  u.Quantity(config['OBJECT']['object_radius_host_star'])
    planet_semi_axis =  u.Quantity(config['OBJECT']['object_semi_major_axis'])
    planet_period = u.Quantity(config['OBJECT']['object_period'])
    planet_inclination = u.Quantity(config['OBJECT']['object_inclination'])
    planet_eccentricity = u.Quantity(config['OBJECT']['object_eccentricity'])
    planet_ephemeris = u.Quantity(config['OBJECT']['object_ephemeris'])
    planet_period = u.Quantity(config['OBJECT']['object_period'])
    planet_omega = u.Quantity(config['OBJECT']['object_omega'])
    ######
    params = batman.TransitParams()       #object to store transit parameters
    params.t0 = 0.                        #time of inferior conjunction
    params.per = 1.                       #orbital period
    # parameter rp #np.sqrt(depth)   #planet radius (in units of stellar radii)
    params.rp = (planet_radius/star_radius).decompose().value # 0.16099487
    params.a = (planet_semi_axis/star_radius).decompose().value # 4.96528796
    params.inc = planet_inclination.to(u.deg).value
    params.ecc = planet_eccentricity.value
    params.w = planet_omega.to(u.deg).value  #longitude of periastron (in degrees)
    params.limb_dark = model_limb_darkening
    params.u = model_limb_darkening_coeff #limb darkening coefficients [u1, u2]
    phase = np.linspace(-model_phase_range/2, model_phase_range/2, model_nphase_points//2)
    #
    return params, phase, planet_period, planet_ephemeris

def define_batman_model(params, phase, atmosphere='atmosphere_transmission.dat'):
    """
    This function reads at table containing the atmosphere transmission versus the wavelength.
    Then it uses the batman package to calculate the
    light curves for each wavelength.
    We assume here a symmetric transit signal, that the
    secondary transit is at phase 0.5 and primary transit at 0.0.
    
    Parameters
    ----------
    params: 'batman parameter object'
        model parameters for batman lightcurve mode (object)
    phase: 'numpy array'
        Orbital phase of planet where the lightcurve is computed
    
    Returns
    -------
    spectra : 'numpy array'
        relative flux function of wavelength and time
    wavelength : 'numpy array'
        wavelength of the different light curves, unit micron
    """
    model = batman.TransitModel(params, phase)
    flux = model.light_curve(params)
    plt.figure()
    plt.plot(phase, flux, '+')
    plt.xlabel("Time from central transit")
    plt.ylabel("Relative flux")
    plt.show()
    #
    print(atmosphere)
    wavelength, transmission = read_depth_ecsv(atmosphere)
    nw = wavelength.size
    nt = phase.size
    #
    plt.figure()
    plt.plot(wavelength, transmission)
    #
    spectra = np.zeros([nw, nt])
    rp = params.rp
    for i in np.arange(nw):
        params.rp = rp*transmission[i]  #updates planet radius
        spectra[i,:] = model.light_curve(params)
    plot_spectrum_2d(spectra, wavelength, phase, 'enfin', normalise=False)
    return spectra, wavelength

def simulate_spectra_correl(file_object='cascade_WASP43b_object.ini', file_model='cascade_WASP43b_calibrate_planet_spectrum.ini', verbose=False):

    #
    #init_path = os.path.join(os.path.dirname(cascade.__path__[0]), 'tests/init_files')
    import inspect
    init_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'init_files')
    print(os.path.join(init_path, file_model))
    config = ConfigObj(os.path.join(init_path, file_model))
    print('OBSERVATIONS', config['OBSERVATIONS'])
    pattern = config['OBSERVATIONS']['observations_id']
    save_path = config['OBSERVATIONS']['observations_path']
    target_name = config['OBSERVATIONS']['observations_target_name']
    pattern = config['OBSERVATIONS']['observations_id']
    #
    spectrum_amplitude = float(config['SIMULATION']['spectrum_amplitude'])
    pick_amplitude = float(config['SIMULATION']['pick_amplitude'])
    atmosphere_file = config['SIMULATION']['atmosphere_file']

    #######  read the ini files and create the batmn object parameters #######
    params, phase, planet_period, planet_ephemeris = return_par_from_ini(file_object=os.path.join(init_path, file_object) ,file_model=os.path.join(init_path, file_model))

    #######  compute the curve ligth  #######
    spectra, wavelength = define_batman_model(params, phase, atmosphere=atmosphere_file)

    lightcurves = spectra*spectrum_amplitude
    noisy_lightcurves = add_noise(lightcurves)
    tso_lightcurves = add_sine_noise(lightcurves, pick_amplitude)
    if(verbose):plot_spectrum_2d(tso_lightcurves, wavelength, phase, 'enfin', normalise=False)

    # data/Generic/Octave/WASP43b/SPECTRA/
    out_dir = os.path.join(save_path, target_name)
    if(config['OBSERVATIONS']['observations_data']=='SPECTRUM'):
        out_dir = os.path.join(out_dir,'SPECTRA')
    #
    uu = u.electron/u.second
    mask = np.full((spectra.shape), False, dtype=bool) # the False are good points
    uncertainties = np.sqrt(lightcurves)
    times = phase*planet_period+planet_ephemeris
    os.makedirs(out_dir,  exist_ok=True)

    write_spectra_dir(out_dir, pattern, wavelength*u.micron, tso_lightcurves*uu, uncertainties*uu,
                          mask, times, verbose=False, overwrite=True)
    return

simulate_spectra_correl(verbose=True)
