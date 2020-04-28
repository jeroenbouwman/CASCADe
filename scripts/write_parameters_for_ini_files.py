#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============
write_parameters_for_ini_files.py
===============

Write the Object Ini Files containing all parameters of the stars and exoplanets needed to run CASCADe

It needs at least 2 catalogs: one with the list of the planets observed by HST (WFC3),
and one with the physical parameters of the stars and exoplanets.

Usage:
python write_parameters_for_ini_files.py

Returns
One object ini file per planet observed by HST

Author:
  Raphael Peralta <raphael.peralta@cea.fr>
Year: 2020
"""

import os
import pandas as pd
import numpy as np
import pickle


def replace_nan_values(values_to_check, values_to_replace):
    """
    Replace missing values in a catalog with the corresponding value in another catalog.
    TODO

    Parameters
    ----------
    values_to_check:
        DESCRIPTION.
    values_to_replace:
        DESCRIPTION.

    Returns
    -------
    new_values :
        DESCRIPTION.

    """
    new_values = np.asarray(values_to_check)
    for i in range(new_values.size):
        if np.isnan(values_to_check[i]):
            new_values[i] = values_to_replace[i]
    return new_values


def remove_space(object_names):
    """
    Remove spaces in planet names.

    Parameters
    ----------
    object_names : list of string of characters
        planet names

    Returns
    -------
    new_name : list of string of characters
        planet names without spaces.

    """
    new_name = np.asarray(object_names)
    for i in range(new_name.size):
        if not pd.isna(object_names[i]):
            new_name[i] = (''.join(new_name[i].split(' ')))
            new_name[i] = (''.join(new_name[i].split('_')))
    return new_name


def remove_binary_name(object_names):
    """
    Remove letter corresponding to binary system in planet names.

    Parameters
    ----------
    object_names : list of string of characters
        planet names

    Returns
    -------
    new_name : list of string of characters
        planet names without binary letters
    """
    new_name = np.asarray(object_names)
    for i in range(new_name.size):
        if not pd.isna(new_name[i]):
            new_name[i] = new_name[i].replace(' (AB) ', '')
            new_name[i] = new_name[i].replace(' A ', '')
            new_name[i] = new_name[i].replace(' B ', '')
    return new_name


def remove_duplicate(object_names):
    """
    Remove duplicates.

    Parameters
    ----------
    object_names : list of string of characters
        planet names

    Returns
    -------
    new_name : list of string of characters
        planet names without duplicates.
    """
    new_name = np.asarray(object_names)
    new_name = np.unique(new_name)
    return new_name


def get_alternative_planet_name(list_planet, list_star):
    """
    add the letter of the planet to one of the names of the host star

    Parameters
    ----------
    list_planet : list of string of characters
        planet name
    list_star : list of string of characters
        star name

    Returns
    -------
    new_planet_name : list of string of characters
        alternative planet name

    """
    new_planet_name = np.asarray(list_star)
    for i in range(new_planet_name.size):
        if not pd.isna(new_planet_name[i]):
            planet_letter = list_planet[i][-1]
            new_planet_name[i] = list_star[i] + planet_letter
    return new_planet_name


def print_parameter_missing(planet_name, planet_radius, star_radius, star_temperature, semi_major_axis, inclination,
                            eccentricity, omega, ephemeris, period):
    """
    print the missing parameters when writing the ini files

    Parameters
    ----------
    planet_name : string
    planet_radius : float
    star_radius : float
    star_temperature : float
    semi_major_axis : float
    inclination : float
    eccentricity : float
    omega : float
    ephemeris : float
    period : float

    Returns
    -------
    print on the shell a Warning message with the missing parameter(s)

    """
    names = np.asarray(['planet radius', 'star radius', 'star temperature', 'semi major axis',
                        'inclination', 'eccentricity', 'omega', 'ephemeris', 'period'])
    values = np.asarray([planet_radius, star_radius, star_temperature, semi_major_axis, inclination,
                         eccentricity, omega, ephemeris, period])
    where_nan = np.isnan(values)
    if names[where_nan].size == 1:
        print("Warning: for the planet {}, the parameter {} is missing".format(planet_name, names[where_nan]))
    if names[where_nan].size > 1:
        print("Warning: for the planet {}, the parameters {} are missing".format(planet_name, names[where_nan]))


def get_catalogue(path_catalogue, name_catalogue, url_catalogue, update=False):
    """
    Read file containing the HST observations
    and extract the list of planet, type of observations and the filter used

    Parameters
    ----------
    path_catalogue : string
        catalogue path
    name_catalogue : string
        catalogue filename. The format have to be a csv file
    url_catalogue  : string
        link to the online catalogue
    update : boolean
        if True, the online catalogue is used

    Returns
    -------
    catalogue_data : pandas object
         full catalogue (last update if update==True)

    """
    if update:
        catalogue_data = pd.read_csv(url_catalogue)
    else:
        catalogue_data = pd.read_csv(os.path.join(path_catalogue, name_catalogue))

    return catalogue_data


def get_hst_parameters(path_catalogue, name_catalogue):
    """
    Read the file containing the HST observations data
    and extract the list of planet, type of observations and the filter used

    Parameters
    ----------
    path_catalogue : string
    name_catalogue : string

    Returns
    -------
    hst_data : pandas object
         planet name, type of observation and filter used

    """
    hst_observations = pickle.load(open(os.path.join(path_catalogue, name_catalogue), 'rb'))

    hst_data = pd.DataFrame({'PLANET': [i['planet'] for i in hst_observations.values()],
                             'OBSERVATION': [i['observation'] for i in hst_observations.values()],
                             'FILTER': [i['filter'] for i in hst_observations.values()]})

    hst_data['PLANET'] = remove_binary_name(hst_data['PLANET'])
    hst_data['PLANET'] = remove_space(hst_data['PLANET'])

    hst_data = hst_data.sort_values(by=['PLANET'])

    return hst_data


def reformat_exoplanets_a_data(exo_a_cat):
    """
    Keeps only the relevant parameters from the *Exoplanet.a* catalogue
    Planet name and alternative planet name are modified
    The header is also modified in order to be uniform with the other catalogue

    Parameters
    ----------
    exo_a_cat : pandas object
        raw catalogue

    Returns
    -------
    exo_a_data : pandas object
         relevant parameters with an uniform header

    """
    exo_a_cat = exo_a_cat.rename({'#NamePlanet': 'PLANET', 'NamePlanet': 'PLANET',
                'RAdeg_Simbad': 'RA', 'DECdeg_SIMBAD': 'DEC', 'RAdeg': 'RA', 'DECdeg': 'DEC',
                'RadiusPlanet_EU': 'RPLANET', 'radius_error_max_EU': 'RPLANETUPPER', 'radius_error_min_EU': 'RPLANETLOWER',
                'semi_major_axis_EU': 'A', 'semi_major_axis_error_max_EU': 'AUPPER', 'semi_major_axis_error_min_EU': 'ALOWER',
                'inclination_EU': 'I', 'inclination_error_max_EU': 'IUPPER', 'inclination_error_min_EU': 'ILOWER',
                'eccentricity_EU': 'ECC', 'eccentricity_error_max_EU': 'ECCUPPER', 'eccentricity_error_min_EU': 'ECCLOWER',
                'omega_EU': 'OM', 'omega_error_max_EU': 'OMUPPER', 'omega_error_min_EU': 'OMLOWER',
                'orbital_period_EU': 'PER', 'orbital_period_error_max_EU': 'PERUPPER', 'orbital_period_error_min_EU': 'PERLOWER',
                'magK_EU': 'KMAG',
                'tzero_tr_EU': 'TT', 'tzero_tr_error_max_EU': 'TTUPPER', 'tzero_tr_error_min_EU': 'TTLOWER',
                'star_name': 'STAR', 'AltName': 'OTHERNAME',
                'star_Radius_EU': 'RSTAR', 'star_eRadius_max_EU': 'RSTARUPPER', 'star_eRadius_min_EU': 'RSTARLOWER',
                'star_Teff_EU': 'TEFF', 'star_eTeff_max_EU': 'TEFFUPPER', 'star_eTeff_min_EU': 'TEFFLOWER',
                'star_Metallicity_EU': 'FE', 'star_eMetallicity_max_EU': 'FEUPPER', 'star_eMetallicity_min_EU': 'FELOWER',
                'loggstar_Ex': 'LOGG'
                }, axis='columns')

    exo_a_cat = exo_a_cat.assign(KMAGUPPER=np.nan, KMAGLOWER=np.nan, LOGGUPPER=np.nan, LOGGLOWER=np.nan)

    exo_a_cat['PLANET'] = remove_binary_name(exo_a_cat['PLANET'])
    exo_a_cat['PLANET'] = remove_space(exo_a_cat['PLANET'])

    exo_a_cat['OTHERNAME'] = remove_binary_name(exo_a_cat['OTHERNAME'])
    exo_a_cat['OTHERNAME'] = remove_space(exo_a_cat['OTHERNAME'])
    exo_a_cat['OTHERNAME'] = get_alternative_planet_name(exo_a_cat['PLANET'], exo_a_cat['OTHERNAME'])

    exo_a_data = exo_a_cat[['PLANET', 'OTHERNAME', 'RA', 'DEC', 'RPLANET', 'RPLANETUPPER', 'RPLANETLOWER',
                            'A', 'ALOWER', 'AUPPER', 'I', 'IUPPER', 'ILOWER', 'ECC', 'ECCUPPER', 'ECCLOWER',
                            'OM', 'OMUPPER', 'OMLOWER', 'PER', 'PERUPPER', 'PERLOWER',
                            'KMAG', 'KMAGUPPER', 'KMAGLOWER', 'TT', 'TTUPPER', 'TTLOWER',
                            'STAR', 'RSTAR', 'RSTARUPPER', 'RSTARLOWER', 'TEFF', 'TEFFUPPER', 'TEFFLOWER',
                            'FE', 'FEUPPER', 'FELOWER', 'LOGG', 'LOGGUPPER', 'LOGGLOWER']]

    exo_a_data = exo_a_data.sort_values(by=['PLANET'])

    return exo_a_data


def reformat_exoplanets_org_data(exo_org_cat):
    """
    Keeps only the relevant parameters from the *Exoplanets.org* catalogue
    Planet name and alternative planet name are modified
    The header is also modified in order to be uniform with the other catalogue

    Parameters
    ----------
    exo_org_cat : pandas object
        raw catalogue

    Returns
    -------
    exo_org_data : pandas object
         relevant parameters with an uniform header

    """
    exo_org_cat = exo_org_cat.rename({'NAME': 'PLANET', 'R': 'RPLANET',
                                      'RUPPER': 'RPLANETUPPER', 'RLOWER': 'RPLANETLOWER'
                                      }, axis='columns')

    exo_org_cat = exo_org_cat.assign(KMAG=np.nan, KMAGUPPER=np.nan, KMAGLOWER=np.nan)

    exo_org_cat['PLANET'] = remove_binary_name(exo_org_cat['PLANET'])
    exo_org_cat['PLANET'] = remove_space(exo_org_cat['PLANET'])

    exo_org_cat['OTHERNAME'] = remove_binary_name(exo_org_cat['OTHERNAME'])
    exo_org_cat['OTHERNAME'] = remove_space(exo_org_cat['OTHERNAME'])
    exo_org_cat['OTHERNAME'] = get_alternative_planet_name(exo_org_cat['PLANET'], exo_org_cat['OTHERNAME'])

    exo_org_data = exo_org_cat[['PLANET', 'OTHERNAME', 'RA', 'DEC', 'RPLANET', 'RPLANETUPPER', 'RPLANETLOWER',
                                'A', 'ALOWER', 'AUPPER', 'I', 'IUPPER', 'ILOWER', 'ECC', 'ECCUPPER', 'ECCLOWER',
                                'OM', 'OMUPPER', 'OMLOWER', 'PER', 'PERUPPER', 'PERLOWER',
                                'KMAG', 'KMAGUPPER', 'KMAGLOWER', 'TT', 'TTUPPER', 'TTLOWER',
                                'STAR', 'RSTAR', 'RSTARUPPER', 'RSTARLOWER', 'TEFF', 'TEFFUPPER', 'TEFFLOWER',
                                'FE', 'FEUPPER', 'FELOWER', 'LOGG', 'LOGGUPPER', 'LOGGLOWER']]

    exo_org_data = exo_org_data.sort_values(by=['PLANET'])

    return exo_org_data


def reformat_nasa_exoplanet_archive_data(nea_cat):
    """
    Keeps only the relevant parameters from the *Nasa Exoplanet Archive* catalogue
    Planet name and alternative planet name are modified
    The header is also modified in order to be uniform with the other catalogue

    Parameters
    ----------
    nea_cat : pandas object
        raw catalogue

    Returns
    -------
    nea_data : pandas object
         relevant parameters with an uniform header

    """
    nea_cat = nea_cat.rename({'pl_name': 'PLANET', 'ra': 'RA', 'dec': 'DEC',
                              'pl_radj': 'RPLANET', 'pl_radjerr1': 'RPLANETUPPER', 'pl_radjerr2': 'RPLANETLOWER',
                              'pl_orbsmax': 'A', 'pl_orbsmaxerr1': 'AUPPER', 'pl_orbsmaxerr2': 'ALOWER',
                              'pl_orbincl': 'I', 'pl_orbinclerr1': 'IUPPER', 'pl_orbinclerr2': 'ILOWER',
                              'pl_orbeccen': 'ECC', 'pl_orbeccenerr1': 'ECCUPPER', 'pl_orbeccenerr2': 'ECCLOWER',
                              'pl_orblper': 'OM', 'pl_orblpererr1': 'OMUPPER', 'pl_orblpererr2': 'OMLOWER',
                              'pl_orbper': 'PER', 'pl_orbpererr1': 'PERUPPER', 'pl_orbpererr2': 'PERLOWER',
                              'st_k': 'KMAG', 'st_kerr': 'KMAGUPPER',
                              'pl_tranmid': 'TT', 'pl_tranmiderr1': 'TTUPPER', 'pl_tranmiderr2': 'TTLOWER',
                              'pl_hostname': 'STAR', 'hd_name': 'OTHERNAME',
                              'st_rad': 'RSTAR', 'st_raderr1': 'RSTARUPPER', 'st_raderr2': 'RSTARLOWER',
                              'st_teff': 'TEFF', 'st_tefferr1': 'TEFFUPPER', 'st_tefferr2': 'TEFFLOWER',
                              'st_metfe': 'FE', 'st_metfeerr1': 'FEUPPER', 'st_metfeerr2': 'FELOWER',
                              'st_logg': 'LOGG', 'st_loggerr1': 'LOGGUPPER', 'st_loggerr2': 'LOGGLOWER'
                              }, axis='columns')

    nea_cat = nea_cat.assign(KMAGLOWER=nea_cat['KMAGUPPER'])

    nea_cat['PLANET'] = remove_binary_name(nea_cat['PLANET'])
    nea_cat['PLANET'] = remove_space(nea_cat['PLANET'])

    nea_cat['OTHERNAME'] = remove_binary_name(nea_cat['OTHERNAME'])
    nea_cat['OTHERNAME'] = remove_space(nea_cat['OTHERNAME'])
    nea_cat['OTHERNAME'] = get_alternative_planet_name(nea_cat['PLANET'], nea_cat['OTHERNAME'])

    nea_data = nea_cat[['PLANET', 'OTHERNAME', 'RA', 'DEC', 'RPLANET', 'RPLANETUPPER', 'RPLANETLOWER',
                        'A', 'ALOWER', 'AUPPER', 'I', 'IUPPER', 'ILOWER', 'ECC', 'ECCUPPER', 'ECCLOWER',
                        'OM', 'OMUPPER', 'OMLOWER', 'PER', 'PERUPPER', 'PERLOWER', 'KMAG', 'KMAGUPPER', 'KMAGLOWER',
                        'TT', 'TTUPPER', 'TTLOWER', 'STAR', 'RSTAR', 'RSTARUPPER', 'RSTARLOWER',
                        'TEFF', 'TEFFUPPER', 'TEFFLOWER', 'FE', 'FEUPPER', 'FELOWER', 'LOGG', 'LOGGUPPER', 'LOGGLOWER']]

    nea_data = nea_data.sort_values(by=['PLANET'])

    return nea_data


def reformat_tepcat_data(tepcat_allplanets_cat, tepcat_observable_cat):
    """
    Keeps only the relevant parameters from the *Tepcat* catalogue
    Planet name and alternative planet name are modified
    The header is also modified in order to be uniform with the other catalogue

    Parameters
    ----------
    tepcat_allplanets_cat : pandas object
        raw catalogue

    Returns
    -------
    tepcat_data : pandas object
         relevant parameters with an uniform header

    """
    tepcat_allplanets_cat[tepcat_allplanets_cat.columns[0]] = \
        remove_binary_name(tepcat_allplanets_cat[tepcat_allplanets_cat.columns[0]])
    tepcat_allplanets_cat[tepcat_allplanets_cat.columns[0]] = \
        remove_space(tepcat_allplanets_cat[tepcat_allplanets_cat.columns[0]])

    tepcat_observable_cat[tepcat_observable_cat.columns[0]] = \
        remove_binary_name(tepcat_observable_cat[tepcat_observable_cat.columns[0]])
    tepcat_observable_cat[tepcat_observable_cat.columns[0]] = \
        remove_space(tepcat_observable_cat[tepcat_observable_cat.columns[0]])

    tepcat_merged_cat = pd.concat([tepcat_allplanets_cat, tepcat_observable_cat], axis=1, sort=True, join='inner')

    tepcat_data = pd.DataFrame({'PLANET': tepcat_merged_cat[tepcat_merged_cat.columns[0]],
                                'OTHERNAME': tepcat_merged_cat[tepcat_merged_cat.columns[43]],
                                'RA': tepcat_merged_cat[tepcat_merged_cat.columns[47]],
                                'DEC': tepcat_merged_cat[tepcat_merged_cat.columns[50]],
                                'RPLANET': tepcat_merged_cat[tepcat_merged_cat.columns[29]],
                                'RPLANETUPPER': tepcat_merged_cat[tepcat_merged_cat.columns[30]],
                                'RPLANETLOWER': tepcat_merged_cat[tepcat_merged_cat.columns[31]],
                                'A': tepcat_merged_cat[tepcat_merged_cat.columns[23]],
                                'AUPPER': tepcat_merged_cat[tepcat_merged_cat.columns[24]],
                                'ALOWER': tepcat_merged_cat[tepcat_merged_cat.columns[25]],
                                'I': '',
                                'IUPPER': '',
                                'ILOWER': '',
                                'ECC': tepcat_merged_cat[tepcat_merged_cat.columns[20]],
                                'ECCUPPER': tepcat_merged_cat[tepcat_merged_cat.columns[21]],
                                'ECCLOWER': tepcat_merged_cat[tepcat_merged_cat.columns[22]],
                                'OM': '',
                                'OMUPPER': '',
                                'OMLOWER': '',
                                'PER': tepcat_merged_cat[tepcat_merged_cat.columns[57]],
                                'PERUPPER': tepcat_merged_cat[tepcat_merged_cat.columns[58]],
                                'PERLOWER': tepcat_merged_cat[tepcat_merged_cat.columns[58]],
                                'KMAG': tepcat_merged_cat[tepcat_merged_cat.columns[52]],
                                'KMAGUPPER': '',
                                'KMAGLOWER': '',
                                'TT': tepcat_merged_cat[tepcat_merged_cat.columns[55]],
                                'TTUPPER': tepcat_merged_cat[tepcat_merged_cat.columns[56]],
                                'TTLOWER': tepcat_merged_cat[tepcat_merged_cat.columns[56]],
                                'STAR': '',
                                'RSTAR': tepcat_merged_cat[tepcat_merged_cat.columns[10]],
                                'RSTARUPPER': tepcat_merged_cat[tepcat_merged_cat.columns[11]],
                                'RSTARLOWER': tepcat_merged_cat[tepcat_merged_cat.columns[12]],
                                'TEFF': tepcat_merged_cat[tepcat_merged_cat.columns[1]],
                                'TEFFUPPER': tepcat_merged_cat[tepcat_merged_cat.columns[2]],
                                'TEFFLOWER': tepcat_merged_cat[tepcat_merged_cat.columns[3]],
                                'FE': tepcat_merged_cat[tepcat_merged_cat.columns[4]],
                                'FEUPPER': tepcat_merged_cat[tepcat_merged_cat.columns[5]],
                                'FELOWER': tepcat_merged_cat[tepcat_merged_cat.columns[6]],
                                'LOGG': tepcat_merged_cat[tepcat_merged_cat.columns[13]],
                                'LOGGUPPER': tepcat_merged_cat[tepcat_merged_cat.columns[14]],
                                'LOGGLOWER': tepcat_merged_cat[tepcat_merged_cat.columns[15]]
                                })

    tepcat_data = tepcat_data.replace(-1, np.nan)

    tepcat_data = tepcat_data.sort_values(by=['PLANET'])

    return tepcat_data


def merge_catalogue(catalogue_a, catalogue_b, verbose=False):
    """
    It returns the merged catalogue between catalogue_a and catalogue_b.
    The merged catalogue contains the planets from catalogue_a with the data from catalogue_b
    (or NaN if no correspondence).
    Both columns 'PLANET' and 'OTHERNAME' from catalogue_b are used for the correspondence with catalogue_a.

    Parameters
    ----------
    catalogue_a : pandas object
        catalogue containing the HST planet names to compare with catalogue_b
    catalogue_b : pandas object
        catalogue containing the planet data
    verbose : boolean
        if True, raise a warning writing the name of the planets without correspondence between both catalogues.

    Returns
    -------
    b_merge_a : pandas object
         new catalogue: merging of catalogue_a and catalogue_b.

    """
    # b_match_a = catalogue_b.loc[catalogue_b['PLANET'].isin(catalogue_a) |
    # (catalogue_b['OTHERNAME'].isin(catalogue_a) & ~catalogue_b['PLANET'].isin(catalogue_a))]

    only_in_othername = catalogue_b['OTHERNAME'].isin(catalogue_a) \
                        & ~catalogue_b['PLANET'].isin(catalogue_a)

    catalogue_b['PLANET'][only_in_othername] = catalogue_b['OTHERNAME'][only_in_othername]

    not_match = catalogue_a.loc[~catalogue_a.isin(catalogue_b['PLANET'])]
    if verbose and not_match.size > 0:
        print('Warning: list of planets with no HST observations: {}'.format(not_match.values))

    b_merge_a = pd.merge(catalogue_a, catalogue_b, how='left', on=['PLANET'])

    return b_merge_a


def check_folder(folder_path, create_directory, verbose=True):
    """
    Check if a folder exists.
    If not, a warning appears (if verbose=True) and the folder is created (if create_directory=True)

    Parameters
    ----------
    folder_path : string
        path and name of the folder
    create_directory : boolean
        if True, the missing folder is created
    verbose : boolean
        if True, raise warnings if the folder do not exist and if it is created.

    Returns
    -------
    None

    """
    if not(os.path.isdir(folder_path)):
        if verbose:
            print("Warning: the directory '{}' doesn't exist.".format(folder_path))
        if create_directory:
            if verbose:
                print("Creating the folder '{}'".format(folder_path))
            os.makedirs(folder_path)


def write_diff(exoplanet_a, exoplanets_org, nea_catalog, tepcat):

    parameter = ['I', 'RPLANET', 'PER', 'TT', 'A']

    writer = pd.ExcelWriter('diff.xlsx', engine='xlsxwriter')

    exoplanet_a.to_excel(writer, sheet_name='exo_a_merged_hst', index=False)
    exoplanets_org.to_excel(writer, sheet_name='exo_org_merged_hst', index=False)
    nea_catalog.to_excel(writer, sheet_name='nea_merged_hst', index=False)
    tepcat.to_excel(writer, sheet_name='tepcat_merged_hst', index=False)

    for param in parameter:

        ratio_a = exoplanet_a[param]
        ratio_org = exoplanets_org[param]
        ratio_nea = nea_catalog[param]
        ratio_tep = tepcat[param]

        if param == 'RPLANET':
            ratio_a = exoplanet_a['RPLANET'] / exoplanet_a['RSTAR']
            ratio_org = exoplanets_org['RPLANET'] / exoplanets_org['RSTAR']
            ratio_nea = nea_catalog['RPLANET'] / nea_catalog['RSTAR']
            ratio_tep = tepcat['RPLANET'] / tepcat['RSTAR']

        if param == 'I':
            ratio_tep = np.full(ratio_tep.size, np.nan)

        org_diff_a = (ratio_org - ratio_a) / ratio_a*100
        nea_diff_a = (ratio_nea - ratio_a) / ratio_a*100
        tep_diff_a = (ratio_tep - ratio_a) / ratio_a*100

        nea_diff_org = (ratio_nea - ratio_org) / ratio_org*100
        tep_diff_org = (ratio_tep - ratio_org) / ratio_org*100

        tep_diff_nea = (ratio_tep - ratio_nea) / ratio_nea*100

        if param == 'TT':
            org_diff_a = (ratio_org - ratio_a) / exoplanet_a['PER']
            nea_diff_a = (ratio_nea - ratio_a) / exoplanet_a['PER']
            tep_diff_a = (ratio_tep - ratio_a) / exoplanet_a['PER']

            nea_diff_org = (ratio_nea - ratio_org) / exoplanets_org['PER']
            tep_diff_org = (ratio_tep - ratio_org) / exoplanets_org['PER']

            tep_diff_nea = (ratio_tep - ratio_nea) / nea_catalog['PER']

        diff = pd.DataFrame({'(ORG-A)/A': org_diff_a, '(NEA-A)/A': nea_diff_a, '(TEP-A)/A': tep_diff_a,
                             '(NEA-ORG/ORG)': nea_diff_org, '(TEP-ORG)/ORG': tep_diff_org,
                             '(TEP-NEA/NEA)': tep_diff_nea})

        diff.to_excel(writer, sheet_name='diff'+param, index=False)

    writer.save()


# Initialisation
verbose = True

#  planets observed by HST WFC3
list_planets_hst_wfc3_path = '../data/archive_databases/HST/WFC3'
list_planets_hst_wfc3_filename = 'WFC3_files.pickle'
object_ini_files_path = 'object_init_files'

# Exoplanet.a catalogue
path_catalogue_exoplanets_a = '../data/exoplanet_data/EXOPLANETS_A'
name_catalogue_exoplanets_a = 'exoplanets_a_full.csv'
url_catalogue_exoplanets_a = 'http://svo2.cab.inta-csic.es/vocats/v2/exostars/cs.php?RA=180.000000&DEC=0.000000&S' \
                             'R=180.000000&VERB=1&objlist=-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,' \
                             '22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,' \
                             '51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,' \
                             '80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,' \
                             '106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,' \
                             '127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,' \
                             '148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,' \
                             '169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,' \
                             '190,191,192,193,194,195,196,197,198,199,200,201&fldlist=-1,3,4,5,6,7,8,9,12,13,14,15,' \
                             '16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,' \
                             '45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,' \
                             '74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,' \
                             '102,103,104,105,106,107,108,109,110,111,112,113&nocoor=1&format=ascii'

# Exoplanets.org catalogue
path_catalogue_exoplanets_org = '../data/exoplanet_data/exoplanets.org/'
name_catalogue_exoplanets_org = 'exoplanets.csv'
url_catalogue_exoplanets_org = 'http://www.exoplanets.org/csv-files/exoplanets.csv'

# Nasa Exoplanet Archive catalogue
path_catalogue_nasa_exoplanet_archive = '../data/exoplanet_data/NASAEXOPLANETARCHIVE/'
name_catalogue_nasa_exoplanet_archive = 'nasaexoplanetarchive_full.csv'
url_catalogue_nasa_exoplanet_archive = ''

# Tepcat allplanet catalogue
path_allplanet_tepcat = '../data/exoplanet_data/tepcat/'
name_allplanet_tepcat = 'allplanets.csv'
url_allplanet_tepcat = 'http://www.astro.keele.ac.uk/jkt/tepcat/allplanets-csv.csv'

# Tepcat observables catalogue
path_observable_tepcat = '../data/exoplanet_data/tepcat/'
name_observable_tepcat = 'observables.csv'
url_observable_tepcat = 'http://www.astro.keele.ac.uk/jkt/tepcat/observables.csv'


# catalogue reading
hst_data = get_hst_parameters(list_planets_hst_wfc3_path, list_planets_hst_wfc3_filename)

exoplanet_a_catalogue = get_catalogue(path_catalogue_exoplanets_a,
                                      name_catalogue_exoplanets_a,
                                      url_catalogue=url_catalogue_exoplanets_a, update=False)

exoplanet_org_catalogue = get_catalogue(path_catalogue_exoplanets_org,
                                        name_catalogue_exoplanets_org,
                                        url_catalogue=url_catalogue_exoplanets_org, update=False)

nea_catalogue = get_catalogue(path_catalogue_nasa_exoplanet_archive,
                              name_catalogue_nasa_exoplanet_archive,
                              url_catalogue=url_catalogue_nasa_exoplanet_archive, update=False)

tepcat_allplanets_catalogue = get_catalogue(path_allplanet_tepcat,
                                            name_allplanet_tepcat,
                                            url_catalogue=url_allplanet_tepcat, update=False)

tepcat_observable_catalogue = get_catalogue(path_observable_tepcat,
                                            name_observable_tepcat,
                                            url_catalogue=url_observable_tepcat, update=False)


# data reading
list_planets_hst_wfc3 = hst_data.drop_duplicates('PLANET')['PLANET']  # remove hst data duplicates

exo_a_merged_hst = reformat_exoplanets_a_data(exoplanet_a_catalogue)
exoplanet_org_data = reformat_exoplanets_org_data(exoplanet_org_catalogue)
nea_data = reformat_nasa_exoplanet_archive_data(nea_catalogue)
tepcat_data = reformat_tepcat_data(tepcat_allplanets_catalogue, tepcat_observable_catalogue)

exo_a_merged_hst = merge_catalogue(list_planets_hst_wfc3, exo_a_merged_hst, verbose=False)
exo_org_merged_hst = merge_catalogue(list_planets_hst_wfc3, exoplanet_org_data, verbose=False)
nea_merged_hst = merge_catalogue(list_planets_hst_wfc3, nea_data, verbose=False)
tepcat_merged_hst = merge_catalogue(list_planets_hst_wfc3, tepcat_data, verbose=False)

# Write merged catalogue with parameters comparaison
#write_diff(exo_a_merged_hst, exo_org_merged_hst, nea_merged_hst, tepcat_merged_hst)


# planet parameters
planet_name = exo_a_merged_hst['PLANET']  # Name of the planet
alternative_planet_name = exo_a_merged_hst['OTHERNAME']  # Other name of the planet
planet_radius = exo_a_merged_hst['RPLANET']  # Radius of the planet [Rjupiter]
semi_major_axis = exo_a_merged_hst['A']  # Orbit semi major axis [au]
inclination = exo_a_merged_hst['I']  # Orbital inclination [degrees]
eccentricity = exo_a_merged_hst['ECC']  # Orbital eccentricity []
omega = exo_a_merged_hst['OM']  # longitude of the ascendent node [degrees]
orbital_period = exo_a_merged_hst['PER']  # Orbital period [days]
kmag = exo_a_merged_hst['KMAG']  # Magnitude in the K filter [magnitude]
ephemeris = exo_a_merged_hst['TT']  # primary transit [JD]

# star parameters
star_name = exo_a_merged_hst['STAR']  # Name of the star
star_radius = exo_a_merged_hst['RSTAR']  # Radius of the star [Rsun]
star_temperature = exo_a_merged_hst['TEFF']  # Star temperature [K]
star_metallicity = exo_a_merged_hst['FE']  # Star metallicity [dex]
star_logg = exo_a_merged_hst['LOGG']  # Star surface gravity [dex(cm/s2)]


# Write object ini files
check_folder(object_ini_files_path, create_directory=True, verbose=verbose)

for i in range(list_planets_hst_wfc3.size):

    object_ini_files_name = 'cascade_' + planet_name[i] + '_object.ini'

    with open(os.path.join(object_ini_files_path, object_ini_files_name), "w") as output_file:
        output_file.write('[OBJECT]\n')
        output_file.write('object_name = {}\n'.format(planet_name[i]))
        output_file.write('object_radius = {} Rjup\n'.format(planet_radius[i]))
        output_file.write('object_radius_host_star = {} Rsun\n'.format(star_radius[i]))
        output_file.write('object_temperature_host_star = {} K\n'.format(star_temperature[i]))
        output_file.write('object_semi_major_axis = {} AU\n'.format(semi_major_axis[i]))
        output_file.write('object_inclination = {} deg\n'.format(inclination[i]))
        output_file.write('object_eccentricity = {}\n'.format(eccentricity[i]))
        output_file.write('object_omega = {} deg\n'.format(omega[i]))
        output_file.write('object_period = {} d\n'.format(orbital_period[i]))
        output_file.write('object_ephemeris = {} d\n'.format(ephemeris[i]))
        output_file.write('object_kmag = {} Kmag\n'.format(kmag[i]))
        output_file.write('object_metallicity_host_star = {} dex\n'.format(star_metallicity[i]))
        output_file.write('object_logg_host_star = {} dex(cm/s2)\n'.format(star_logg[i]))
        output_file.write('\n')
        output_file.write('[CATALOG]\n')
        output_file.write('catalog_use_catalog = False\n')
        output_file.write('catalog_name = MIXTE\n')
        output_file.write('catalog_update = True\n')

    if verbose:
        print_parameter_missing(planet_name[i], planet_radius[i], star_radius[i], star_temperature[i],
                            semi_major_axis[i], inclination[i], eccentricity[i], omega[i],
                            ephemeris[i], orbital_period[i])