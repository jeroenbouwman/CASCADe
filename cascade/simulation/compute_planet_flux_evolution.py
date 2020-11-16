#!/usr/bin/env python
# -*- coding: utf-8 -*-

# THIS PACKAGE COMPUT THE DEPTH OF THE TRANSIT WITH A SIMPLE GEOMETRIC MODEL
# SEE MARINE MARTIN-LAGARDE EXONOODLE FOR A MORE PRECISE TREATMENT

import numpy as np
import math
from astropy import units as u
import matplotlib.pyplot as plt
from configobj import ConfigObj
from cascade.utilities.readwrite_spectra import plot_spectrum_2d

"""
    Compute the flux of the planet or a given phase from the maximum planet flux night and the maximum planet flux day.
    The flux of a planet, for a given phase, is the sum of these 2 terms
    planet_flux_night = (1+cos(phase))/2 * planet_flux_night_0
    planet_flux_day   = (1-cos(phase))/2 * planet_flux_day_0
    For phase = 0 we have only night , which is normal because it is the middle of the transit
  
    :param quantity planet_flux_day_0: amplitude of the flux when all the planet is lit for the observer, with unit usually electron/second
    
    :param quantity planet_flux_nigth_0: amplitude of the flux when all the planet is in the shade for the observer, with unit usually electron/second
    
    :History:
    09 Nov 2019: Created
    @author: Rene Gastaud (CEA-Saclay IRFU/DEDI)
"""

def flux_planet(planet_flux_day, planet_flux_night, angle):
    nw = planet_flux_day.size
    nt = angle.size
    angle = angle.reshape([1,nt])
    result = planet_flux_day.reshape([nw,1])*(1-np.cos(angle))/2. + \
           planet_flux_night.reshape([nw,1])*(1+np.cos(angle))/2.
    return result

def return_par_from_ini_bis(file_object, file_model, verbose=False):
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
    albedo = float(configM['SIMULATION']['planet_bond_albedo'])
    if(verbose): print("albedo ", albedo)
    #######
    configO = ConfigObj(file_object)
    if(verbose): print('config object', configO['OBJECT'])
    planet_name = configO['OBJECT']['object_name']
    planet_radius = u.Quantity(configO['OBJECT']['object_radius'])
    planet_semi_axis =  u.Quantity(configO['OBJECT']['object_semi_major_axis'])
    #planet_period = u.Quantity(configO['OBJECT']['object_period'])
    return planet_radius, planet_semi_axis, albedo
######

#####################################################################################################

def compute_planet_flux_evolution(star_flux,  planet_flux_day_0, planet_flux_night_0,\
                          phase, file_object, file_model, verbose=False):

    """
    Compute the flux of a planet.
    The emitteed flux of a planet, for a given phase is the sum of these 2 terms
        planet_flux_night = (1+cos(phase))/2 * planet_flux_night_0
        planet_flux_day   = (1-cos(phase))/2 * planet_flux_day_0
    For phase = 0 we have only night , which is normal because it is the middle of the transit
    The reflected flux of a planet is given by
        planet_reflected_flux = star_flux*reflectance*((1-cos(phase))/2
    
    BEWARE : for reflection the albedo  is treated as day emission (Lambert ?))
    
    :param quantity star_flux:  flux of the star with unit usually electron/second

    :param quantity planet_flux_day_0: amplitude of the flux when all the planet is lit for the observer, with unit usually electron/second
    
    :param quantity planet_flux_nigth_0: amplitude of the flux when all the planet is in the shade for the observer, with unit usually electron/second

    :param np.array phase : between -0.5 and 0.5
    
    file_object : 'str'
    name of the object configuration file,
    file_model : 'str'
    name of the model configuration file
    
    :return: planet_flux
    :rtype: np.array


    :History:
    15 May 2020 rewritten again !
    @author: Rene Gastaud (CEA-Saclay IRFU/DEDI)
    """


    ##########
    planet_radius, planet_semi_axis, albedo = return_par_from_ini_bis(file_object, file_model, verbose=False)
    ##########  code
    # see  Giuseppe_Morello_Laurea_Magistrale_Fisica December 2012, formula 1.4
    ratio = (planet_radius/planet_semi_axis/2).decompose()  # version 1.1
    reflectance = albedo*(ratio.value**2)
    planet_flux_day_1  = planet_flux_day_0 + reflectance*star_flux
 
    if (verbose):  print("reflectance={:8.6E}".format(reflectance))
    angle = 2*math.pi*phase # *u.radian
    nt = phase.size
    planet_flux2d = flux_planet(planet_flux_day_1, planet_flux_night_0, angle)
    #if(plot): plot_spectrum_2d(planet_flux2d.value, wavelength, phase, 'planet flux', normalise=False)

    return planet_flux2d


