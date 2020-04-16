#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from typing import Any, Union

import pandas as pd
import numpy as np
import pickle


def replace_nan_values(values_to_check, values_to_replace):
    new_values = np.asarray(values_to_check)
    for i in range(new_values.size):
        if np.isnan(values_to_check[i]):
            new_values[i] = values_to_replace[i]
    return new_values


def remove_space(object_names):
    new_name = np.asarray(object_names)
    for i in range(new_name.size):
        new_name[i] = (''.join(object_names[i].split(' ')))
    return new_name


def remove_binary_name(object_names):
    new_name = np.asarray(object_names)
    for i in range(new_name.size):
        new_name[i] = new_name[i].replace(' (AB) ', '')
        new_name[i] = new_name[i].replace(' A ', '')
    return new_name


def remove_duplicate(object_names):
    new_name = np.asarray(object_names)
    new_name = np.unique(new_name)
    return new_name


def get_alternative_planet_name(list_planet, list_star):
    new_planet_name = np.asarray(list_star)
    for i in range(new_planet_name.size):
        planet_letter = list_planet[i][-1]
        new_planet_name[i] = list_star[i] + planet_letter
    return new_planet_name


def print_parameter_missing(planet_name, planet_radius, star_radius, star_temperature, semi_major_axis, inclination,
                            eccentricity, omega, ephemeris):
    names = np.asarray(['planet radius', 'star radius', 'star temperature', 'semi major axis',
                        'inclination', 'eccentricity', 'omega', 'ephemeris'])
    values = np.asarray([planet_radius, star_radius, star_temperature, semi_major_axis, inclination,
                         eccentricity, omega, ephemeris])
    where_nan = np.isnan(values)
    if names[where_nan].size == 1:
        print("Warning: for the planet %s, the parameter %s is missing" % (planet_name, names[where_nan]))
    if names[where_nan].size > 1:
        print("Warning: for the planet %s, the parameters %s are missing" % (planet_name, names[where_nan]))

# Initialisation
verbose = True
list_planets_hst_wfc3_path = '../data/data/HST/archive_data/'
list_planets_hst_wfc3_filename = 'WFC3_files.pickle'
planet_physical_parameters_path = '../data/exoplanet_data/EXOPLANETS_A/'
planet_physical_parameters_filename = 'exoplanets_a_full.csv'
object_ini_files_path = 'object_init_files/'


# Data reading
hst_data = pickle.load(open(list_planets_hst_wfc3_path+list_planets_hst_wfc3_filename, 'rb'))
planet_physical_parameters = pd.read_csv(planet_physical_parameters_path + planet_physical_parameters_filename)


# list of planets observed by HST WFC3
list_planets_hst_wfc3 = [i['planet'] for i in hst_data.values()]

# planet parameters
planet_name = planet_physical_parameters['NamePlanet']  # Name of the planet
# Radius of the planet [Rjupiter] from either Exoplanet.eu or Exeter Librairy
planet_radius = replace_nan_values(planet_physical_parameters['RadiusPlanet_EU'], planet_physical_parameters['Rpl_Ex'])
semi_major_axis = planet_physical_parameters['semi_major_axis_EU']  # Orbit semi major axis [au]
inclination = planet_physical_parameters['inclination_EU']  # Orbital inclination [degrees]
eccentricity = planet_physical_parameters['eccentricity_EU']  # Orbital eccentricity []
omega = planet_physical_parameters['omega_EU']  # longitude of the ascendent node [degrees]
orbital_period = planet_physical_parameters['orbital_period_EU']  # Orbital period [days]
kmag = planet_physical_parameters['magK_EU']  # Magnitude in the K filter [magnitude]
ephemeris = planet_physical_parameters['tzero_tr_EU']  # primary transit [JD]

# star parameters
star_name = planet_physical_parameters['star_name']  # Name of the star
alternative_star_name = planet_physical_parameters['AltName']  # Alternative Name of the star
# Radius of the star [Rsun] from either Exoplanet.eu or Exeter Librairy
star_radius = replace_nan_values(planet_physical_parameters['star_Radius_EU'], planet_physical_parameters['Rstar_Ex'])
star_temperature = planet_physical_parameters['star_Teff_EU']  # Star temperature [K]
star_metallicity = planet_physical_parameters['star_Metallicity_EU']  # Star metallicity [dex]
star_logg = planet_physical_parameters['loggstar_Ex']  # Star surface gravity [dex(cm/s2)]


# Reformate data
planet_name = remove_binary_name(planet_name)
planet_name = remove_space(planet_name)

list_planets_hst_wfc3 = remove_duplicate(list_planets_hst_wfc3)
list_planets_hst_wfc3 = remove_binary_name(list_planets_hst_wfc3)
list_planets_hst_wfc3 = remove_space(list_planets_hst_wfc3)

alternative_star_name = remove_binary_name(alternative_star_name)
alternative_star_name = remove_space(alternative_star_name)
alternative_planet_name = get_alternative_planet_name(planet_name, alternative_star_name)


# Write object ini files
for planet_HST_WFC3 in list_planets_hst_wfc3:

    for whole_planets in range(planet_name.size):
        match = False
        if planet_HST_WFC3 == planet_name[whole_planets] or planet_HST_WFC3 == alternative_planet_name[whole_planets]:

            object_ini_files_name = 'cascade_' + planet_HST_WFC3 + '_object.ini'
            with open(object_ini_files_path + object_ini_files_name, "w") as output_file:
                output_file.write('[OBJECT]\n')
                output_file.write('object_name = %s\n' % planet_name[whole_planets])
                output_file.write('object_radius = %s Rjup\n' % planet_radius[whole_planets])
                output_file.write('object_radius_host_star = %s Rsun\n' % star_radius[whole_planets])
                output_file.write('object_temperature_host_star = %s K\n' % star_temperature[whole_planets])
                output_file.write('object_semi_major_axis = %s AU\n' % semi_major_axis[whole_planets])
                output_file.write('object_inclination = %s deg\n' % inclination[whole_planets])
                output_file.write('object_eccentricity = %s\n' % eccentricity[whole_planets])
                output_file.write('object_omega = %s deg\n' % omega[whole_planets])
                output_file.write('object_period = %s d\n' % orbital_period[whole_planets])
                output_file.write('object_ephemeris = %s d\n' % ephemeris[whole_planets])
                output_file.write('object_kmag = %s Kmag\n' % kmag[whole_planets])
                output_file.write('object_metallicity_host_star = %s dex\n' % star_metallicity[whole_planets])
                output_file.write('object_logg_host_star = %s dex(cm/s2)\n' % star_logg[whole_planets])
                output_file.write('\n')
                output_file.write('[CATALOG]\n')
                output_file.write('catalog_use_catalog = False\n')
                output_file.write('catalog_name = EXOPLANET.EU\n')
                output_file.write('catalog_update = True\n')

            if verbose:
                print_parameter_missing(planet_name[whole_planets], planet_radius[whole_planets],
                                        star_radius[whole_planets], star_temperature[whole_planets],
                                        semi_major_axis[whole_planets], inclination[whole_planets],
                                        eccentricity[whole_planets], omega[whole_planets], ephemeris[whole_planets])

            match = True
            break

    if verbose and not match:
        print('Warning: Planet %s observed by HST with no parameters' % planet_HST_WFC3)
