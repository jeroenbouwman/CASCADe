#  Compute  Star Flux with Batman

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
from cascade.utilities.readwrite_spectra import write_depth_ecsv
from cascade.utilities.readwrite_spectra import write_spectra_dir
from cascade.utilities.noises import add_noise
#from cascade.utilities.noises import add_sine_noise

def add_pick(depths, amplitude, width, location):
    '''
        Add a triangular pick to your depths to simulate and simple variation. Features of the pick are its amplitude,
        its width and its location in the array of depths.
        
        :param depths: Depths of light curves for different wavelengths.
        :param amplitude: Amplitude of the pick.
        :param width: Width of the pick.
        :param location: Location of the apex of the pick.
        :type depths: ndarray: dim: nw; value: percent
        :type amplitude: float; value: percent
        :type width: int > 0
        :type location: int; 0 <= location < nw
        
        :return: depths + pick
        :rtype: ndarray: dim: nw; value: percent
        '''
    pick = np.zeros(2 * width + 1)
    len_p = pick.size
    idx = np.arange(len_p)
    pick[:width] = (amplitude / width) * idx[:width]
    pick[width:] = amplitude * ((2 * width) / width) - (amplitude / width) * idx[width:]
    depths[location - width:location + width + 1] += pick
    return depths

def get_atmosphere_transmission(wavelength, file_model, plot=False, verbose=False ):
    config = ConfigObj(file_model)
    pick_amplitude = float(config['SIMULATION']['pick_amplitude'])
    pick_width = int(config['SIMULATION']['pick_width'])
    trans_file = config['SIMULATION']['atmosphere_file']
    if(verbose):print("file={}, pick_amplitude={:.2e}, pick_width={}".format(trans_file, pick_amplitude, pick_width))
    nw = wavelength.size
    if (pick_amplitude > 0):
        depths = np.ones(nw) * 1
        depths = add_pick(depths, pick_amplitude, pick_width,nw//3)
        depths = add_pick(depths, -pick_amplitude, pick_width, 2*nw//3)
        write_depth_ecsv(trans_file, depths, wavelength)
    else:
        wavelengthA, absorption = read_depth_ecsv(trans_file)
        #interpolate
    #
    if (plot):
        plt.figure()
        plt.title("exoplanet atmosphere transmission")
        plt.plot(wavelength, depths)
        plt.xlabel("wavelength "+str(wavelength.unit))
        plt.ylabel("transmission")
    return depths


#######################

def return_par_from_ini(file_object='cascade_WASP43b_object.ini',
                        file_model='cascade_WASP43b_calibrate_planet_spectrum.ini', verbose=False):
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
    if(verbose): print('MODEL', config['MODEL'])
    #model_phase_range = float(config['MODEL']['model_phase_range'])
    #model_nphase_points = int(config['MODEL']['model_nphase_points'])
    model_limb_darkening=config['MODEL']['model_limb_darkening']
    model_limb_darkening_coeff=config['MODEL']['model_limb_darkening_coeff']
    chaine = model_limb_darkening_coeff[0]+','+model_limb_darkening_coeff[1]
    model_limb_darkening_coeff = ast.literal_eval(chaine)
    if(verbose): print(model_limb_darkening)
    #######
    config = ConfigObj(file_object)
    if(verbose): print('config object', config['OBJECT'])
    planet_name = config['OBJECT']['object_name']
    planet_radius = u.Quantity(config['OBJECT']['object_radius'])
    star_radius =  u.Quantity(config['OBJECT']['object_radius_host_star'])
    planet_semi_axis =  u.Quantity(config['OBJECT']['object_semi_major_axis'])
    planet_period = u.Quantity(config['OBJECT']['object_period'])
    planet_inclination = u.Quantity(config['OBJECT']['object_inclination'])
    planet_eccentricity = u.Quantity(config['OBJECT']['object_eccentricity'])
    #planet_ephemeris = u.Quantity(config['OBJECT']['object_ephemeris'])
    #planet_period = u.Quantity(config['OBJECT']['object_period'])
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
    #
    return params

def define_batman_model(params, phase, atmosphere_transmission):
    """
    This function uses the square root of atmosphere transmission for the exoplanet relative radius.
    The atmosphere transmission is a 1D array, same size than wavelength.
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
    atmosphere_transmission 'numpy array'
        Relative variation of the surface of the exoplanet (or transit depth)
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
    #print(atmosphere)
    #wavelength, transmission = read_depth_ecsv(atmosphere)
    nw = atmosphere_transmission.size
    nt = phase.size
    ##
    spectra = np.zeros([nw, nt])
    rp = params.rp
    for i in np.arange(nw):
        params.rp = rp*np.sqrt(atmosphere_transmission[i])  #updates planet radius
        spectra[i,:] = model.light_curve(params)
    #plot_spectrum_2d(spectra, wavelength, phase, 'enfin', normalise=False)
    return spectra

def compute_star_flux_evolution(star_flux1d, wavelength, phase, file_object='cascade_WASP43b_object.ini', file_model='cascade_WASP43b_calibrate_planet_spectrum.ini', verbose=False):

    nt = phase.size
    nw = wavelength.size
    
    #######  read the ini files and create the batmn object parameters #######
    params = return_par_from_ini(file_object=file_object, file_model=file_model)
    
    atmosphere_transmission= \
        get_atmosphere_transmission(wavelength, file_model, plot=False )
    
    #######  compute the  light curves  #######
    spectra = define_batman_model(params, phase, atmosphere_transmission)
                    
    star_flux2d = spectra*star_flux1d.reshape([nw, 1])

    return  star_flux2d

#simulate_spectra_correl(verbose=True)
