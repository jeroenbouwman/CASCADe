#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from typing import Any, Union

import pandas as pd
import numpy as np
import math


def replace_nan_values(values_to_check, values_to_replace):
    new_values = np.asarray(values_to_check)
    for i in range(new_values.size):
        if values_to_check[i] == math.nan:
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


def alternative_planet_name(list_planet, list_star):
    new_planet_name = np.asarray(list_star)
    for i in range(new_planet_name.size):
        planet_letter = list_planet[i][-1]
        new_planet_name[i] = list_star[i] + planet_letter
    return new_planet_name


# Initialisation
list_planets_HST_WFC3_path = ''
list_planets_HST_WFC3_filename = 'catalogue_angelos_planets_hst_wfc3.txt'
planet_physical_parameters_path = ''
planet_physical_parameters_filename = 'catalogue_simon_planets_parameters.csv'
object_ini_files_path = 'object_init_files/'


# Data reading
list_planets_HST_WFC3 = np.loadtxt(list_planets_HST_WFC3_path + list_planets_HST_WFC3_filename, dtype=np.str)
planet_physical_parameters = pd.read_csv(planet_physical_parameters_path + planet_physical_parameters_filename)


# planet parameters
planet_name = planet_physical_parameters['NamePlanet']  # Name of the planet
# Radius of the planet [Rjupiter] from either Exoplanet.eu or Exeter Librairy
planet_radius = replace_nan_values(planet_physical_parameters['RadiusPlanet_EU'], planet_physical_parameters['Rpl_Ex'])
semi_major_axis = planet_physical_parameters['semi_major_axis_EU']  # Semi major axis [au]
inclination = planet_physical_parameters['inclination_EU']  # Orbital inclination [degrees]
eccentricity = planet_physical_parameters['eccentricity_EU']  # Orbital eccentricity []
omega = planet_physical_parameters['omega_EU']  # longitude of the ascendent node [degrees]
rotation = planet_physical_parameters['Period_rotation']  # Rotation period of the star [days]
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

list_planets_HST_WFC3 = remove_duplicate(list_planets_HST_WFC3)
list_planets_HST_WFC3 = remove_binary_name(list_planets_HST_WFC3)
list_planets_HST_WFC3 = remove_space(list_planets_HST_WFC3)

alternative_star_name = remove_binary_name(alternative_star_name)
alternative_star_name = remove_space(alternative_star_name)


# Planet observed by HST
# HST_observation = np.where(mission == 'H')[0]


# Write object ini files
for planet_HST_WFC3 in list_planets_HST_WFC3:
    for whole_planets in range(planet_name.size):
        if planet_HST_WFC3 == planet_name[whole_planets]:
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
                output_file.write('object_period = %s d\n' % rotation[whole_planets])
                output_file.write('object_ephemeris = %s d\n' % ephemeris[whole_planets])
                output_file.write('object_kmag = %s Kmag\n' % kmag[whole_planets])
                output_file.write('object_metallicity_host_star = %s dex\n' % star_metallicity[whole_planets])
                output_file.write('object_logg_host_star = %s dex(cm/s2)\n' % star_logg[whole_planets])
                output_file.write('\n')
                output_file.write('[CATALOG]\n')
                output_file.write('catalog_use_catalog = False\n')
                output_file.write('catalog_name = EXOPLANET.EU\n')
                output_file.write('catalog_update = True\n')

            if  np.isnan(planet_radius[whole_planets]) or  np.isnan(star_radius[whole_planets]) or \
                np.isnan(star_temperature[whole_planets]) or np.isnan(semi_major_axis[whole_planets]) or \
                np.isnan(inclination[whole_planets]) or np.isnan(eccentricity[whole_planets]) or \
                np.isnan(omega[whole_planets]) or np.isnan(ephemeris[whole_planets]) :
                    print('Warning: missing some fundamental parameters for the planet %s' % planet_name[whole_planets])

            break

    if whole_planets == len(planet_name)-1:
        print('Warning: planet %s observed by HST with no parameters' % planet_HST_WFC3)
