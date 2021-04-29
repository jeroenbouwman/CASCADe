#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===============
check_phase_error.py
===============

This code computes and compares the error phase between different catalogs.
The error phase informs if a catalog gives better ephemerids and period values. The lower the value of the error phase,
the more likely it is that the value of the ephemeris and the period is better.

The code reads the parameters of the catalogs and writes the relevant parameters in an excel file
(one sheet per catalog). last sheets are dedicated to the comparison of the error phase and ephemerids values.

Usage:
python check_phase_error.py

Returns
One excel file with several sheets

Author:
  Raphael Peralta <raphael.peralta@cea.fr>
Year: 2020
"""

import os
import pandas as pd
import numpy as np
import astropy.units as u


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


def get_catalogue(path_catalogue, name_catalogue, url_catalogue, skiprows=None, update=False):
    """
    Read file containing the HST observations and extract the parameters

    Parameters
    ----------
    path_catalogue : string
        catalogue path
    name_catalogue : string
        catalogue filename. The format have to be a csv file
    url_catalogue : string
        link to the online catalogue
    skiprows : list-like, int or callable, optional
        skip first rows. By default: None
    update : boolean
        if True, the online catalogue is used

    Returns
    -------
    catalogue_data : pandas object
         full catalogue (last update if update==True)

    """
    if update:
        catalogue_data = pd.read_csv(url_catalogue, skiprows=skiprows)
    else:
        catalogue_data = pd.read_csv(os.path.join(path_catalogue, name_catalogue), skiprows=skiprows)

    return catalogue_data


def get_jeroen_parameters(path_catalogue, name_catalogue):
    """
    Read the catalogue obtained with the "build_local_hst_archive.py" script and extract the parameters

    Parameters
    ----------
    path_catalogue : string
    name_catalogue : string

    Returns
    -------
    jeroen_cat : pandas object
    """
    catalog = pd.read_csv(os.path.join(path_catalogue, name_catalogue), sep=';')

    data = catalog.rename({'planet_name': 'PLANET'}, axis='columns')

    data['PLANET'] = remove_binary_name(data['PLANET'])
    data['PLANET'] = remove_space(data['PLANET'])

    jeroen_cat = data.sort_values(by=['PLANET'])

    return jeroen_cat


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
                             'OBSERVATION': [i['observation'] for i in hst_observations.values()]})

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
    nea_cat = nea_cat.rename({'pl_name': 'PLANET', 'mpl_name': 'PLANET',
                              'ra': 'RA', 'dec': 'DEC',
                              'pl_radj': 'RPLANET', 'pl_radjerr1': 'RPLANETUPPER', 'pl_radjerr2': 'RPLANETLOWER',
                              'mpl_radj': 'RPLANET', 'mpl_radjerr1': 'RPLANETUPPER', 'mpl_radjerr2': 'RPLANETLOWER',
                              'pl_orbsmax': 'A', 'pl_orbsmaxerr1': 'AUPPER', 'pl_orbsmaxerr2': 'ALOWER',
                              'mpl_orbsmax': 'A', 'mpl_orbsmaxerr1': 'AUPPER', 'mpl_orbsmaxerr2': 'ALOWER',
                              'pl_orbincl': 'I', 'pl_orbinclerr1': 'IUPPER', 'pl_orbinclerr2': 'ILOWER',
                              'mpl_orbincl': 'I', 'mpl_orbinclerr1': 'IUPPER', 'mpl_orbinclerr2': 'ILOWER',
                              'pl_orbeccen': 'ECC', 'pl_orbeccenerr1': 'ECCUPPER', 'pl_orbeccenerr2': 'ECCLOWER',
                              'mpl_orbeccen': 'ECC', 'mpl_orbeccenerr1': 'ECCUPPER', 'mpl_orbeccenerr2': 'ECCLOWER',
                              'pl_orblper': 'OM', 'pl_orblpererr1': 'OMUPPER', 'pl_orblpererr2': 'OMLOWER',
                              'mpl_orblper': 'OM', 'mpl_orblpererr1': 'OMUPPER', 'mpl_orblpererr2': 'OMLOWER',
                              'pl_orbper': 'PER', 'pl_orbpererr1': 'PERUPPER', 'pl_orbpererr2': 'PERLOWER',
                              'mpl_orbper': 'PER', 'mpl_orbpererr1': 'PERUPPER', 'mpl_orbpererr2': 'PERLOWER',
                              'st_k': 'KMAG', 'st_kerr': 'KMAGUPPER',
                              'sy_kmag':'KMAG', 'sy_kmagerr1': 'KMAGUPPER', 'sy_kmagerr2': 'KMAGLOWER',
                              'pl_tranmid': 'TT', 'pl_tranmiderr1': 'TTUPPER', 'pl_tranmiderr2': 'TTLOWER',
                              'mpl_tranmid': 'TT', 'mpl_tranmiderr1': 'TTUPPER', 'mpl_tranmiderr2': 'TTLOWER',
                              'pl_hostname': 'STAR', 'mpl_hostname': 'STAR', 'hostname': 'STAR', 'hd_name': 'OTHERNAME',
                              'st_rad': 'RSTAR', 'st_raderr1': 'RSTARUPPER', 'st_raderr2': 'RSTARLOWER',
                              'mst_rad': 'RSTAR', 'mst_raderr1': 'RSTARUPPER', 'mst_raderr2': 'RSTARLOWER',
                              'st_teff': 'TEFF', 'st_tefferr1': 'TEFFUPPER', 'st_tefferr2': 'TEFFLOWER',
                              'mst_teff': 'TEFF', 'mst_tefferr1': 'TEFFUPPER', 'mst_tefferr2': 'TEFFLOWER',
                              'st_metfe': 'FE', 'st_metfeerr1': 'FEUPPER', 'st_metfeerr2': 'FELOWER',
                              'st_met': 'FE', 'st_meterr1': 'FEUPPER', 'st_meterr2': 'FELOWER',
                              'mst_metfe': 'FE', 'mst_metfeerr1': 'FEUPPER', 'mst_metfeerr2': 'FELOWER',
                              'st_logg': 'LOGG', 'st_loggerr1': 'LOGGUPPER', 'st_loggerr2': 'LOGGLOWER',
                              'mst_logg': 'LOGG', 'mst_loggerr1': 'LOGGUPPER', 'mst_loggerr2': 'LOGGLOWER',
                              'ttv_flag': 'TTV'
                              }, axis='columns')

    if 'KMAGLOWER' not in nea_cat:
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
                        'TEFF', 'TEFFUPPER', 'TEFFLOWER', 'FE', 'FEUPPER', 'FELOWER',
                        'LOGG', 'LOGGUPPER', 'LOGGLOWER', 'TTV']]

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
    It returns the merged catalogue between the catalogue_a and catalogue_b.
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

    catalogue_b.loc[only_in_othername, 'PLANET'] = catalogue_b['OTHERNAME'][only_in_othername]

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
            os.mkdir(folder_path)


def write_diff(exoplanet_a, exoplanets_org, nea_catalog, tepcat, output_filename):
    """
    Write the comparison of values of 5 parameters ('I', 'RPLANET', 'PER', 'TT', 'A') between different catalog

    Parameters
    ---------
    exoplanet_a: pandas object
    exoplanets_org: pandas object
    nea_catalog: pandas object
    tepcat: pandas object

    Returns
    -------
    output_filename: Excel file
    """
    parameter = ['I', 'RPLANET', 'PER', 'TT', 'A']

    writer = pd.ExcelWriter(output_filename + '.xlsx', engine='xlsxwriter')

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


def compute_number_period(period, ephemeris, expstart, exptime, verbose=True):
    """
    Computes the number of period

    Parameters
    ----------
    period: float
    ephemeris: float
    expstart: float
    exptime: float
    verbose: boolean

    Returns
    ----------
    nbr_period: float
        number of period
    """
    if verbose:
        if (exptime < period):
            print('Good :observation covers less than one period')
        else:
            print('Bad : observation covers more than one period', exptime/period)

    delta_t = expstart + 2400000.5*u.day - ephemeris    # 'MJD = JD - 2400000.5'
    nbr_period = (delta_t / period).decompose().value

    if verbose:
        print('delta_t = ', delta_t, 'number of periods = ', np.round(nbr_period))

    return np.round(nbr_period)


def compute_phase_error(period, period_error, nbr_period, verbose=False):
    """
    Computes the phase error

    Parameters
    ----------
    period: float
    period_error: float
    nbr_period: float
    verbose: boolean

    Returns
    -------
    phase_error: float

    """
    phase_error = (period_error * nbr_period / period).decompose().value

    if verbose:
        print('phase error = ', phase_error)

    return phase_error


def get_phase_error(visit_cat, param_cat, best_value_only=True):
    """
    Computes the phase error

    Parameters
    ----------
    visit_cat: pandas object
    param_cat: pandas object
    best_value_only: boolean

    Returns
    -------
    results: pandas object
    """
    visit_planet = visit_cat['PLANET'].values
    visit = visit_cat['visit'].values
    observation = visit_cat['observation'].values
    expstart = visit_cat['expstart'].values * u.day
    exptime = visit_cat['exptime'].values * u.s

    param_planet = param_cat['PLANET'].values
    ephemeris = param_cat['TT'].values * u.day  # primary transit [JD]
    ephemeris_err = param_cat['TTUPPER'].values * u.day  # primary transit [JD]
    period = param_cat['PER'].values * u.day
    period_err = param_cat['PERUPPER'].values * u.day

    result = []

    for i in range(visit.size):
        j = np.where(visit_planet[i] == param_planet)[0]

        nbr_period = compute_number_period(period[j], ephemeris[j], expstart[i], exptime[i], verbose=False)
        phase_error = compute_phase_error(period[j], period_err[j], nbr_period, verbose=False)

        if best_value_only:
            try:
                k = np.nanargmin(abs(phase_error))  # nanargmin ignores Nan
            except:
                k = 0   # if phase_error contains only Nan values, returns the first item (a NaN value)

            if 'TTV' in param_cat:
                result.append({'planet_name': visit_planet[i], 'visit': visit[i], 'observation': observation[i],
                           'expstart': expstart[i].value, 'exptime': exptime[i].value,
                           'ephemeris': ephemeris[j].value[k], 'ephemeris_error': ephemeris_err[j].value[k],
                           'period': period[j].value[k], 'period_error': period_err[j].value[k],
                           'nbr_period': nbr_period[k], 'phase_error': phase_error[k],
                           'ttv': param_cat['TTV'][j].values[k]
                           })

            else:
                result.append({'planet_name': visit_planet[i], 'visit': visit[i], 'observation': observation[i],
                           'expstart': expstart[i].value, 'exptime': exptime[i].value,
                           'ephemeris': ephemeris[j].value[k], 'ephemeris_error': ephemeris_err[j].value[k],
                           'period': period[j].value[k], 'period_error': period_err[j].value[k],
                           'nbr_period': nbr_period[k], 'phase_error': phase_error[k]
                           })

        else:

            if 'TTV' in param_cat:
                result.append({'planet_name': visit_planet[i], 'visit': visit[i], 'observation': observation[i],
                           'expstart': expstart[i].value, 'exptime': exptime[i].value,
                           'ephemeris': ephemeris[j].value, 'ephemeris_error': ephemeris_err[j].value,
                           'period': period[j].value, 'period_error': period_err[j].value,
                           'nbr_period': nbr_period, 'phase_error': phase_error,
                           'ttv': param_cat['TTV'][j].values
                           })

            else:
                result.append({'planet_name': visit_planet[i], 'visit': visit[i], 'observation': observation[i],
                           'expstart': expstart[i].value, 'exptime': exptime[i].value,
                           'ephemeris': ephemeris[j].value, 'ephemeris_error': ephemeris_err[j].value,
                           'period': period[j].value, 'period_error': period_err[j].value,
                           'nbr_period': nbr_period, 'phase_error': phase_error
                           })

    result = pd.DataFrame(result)

    return result


def main(root_dir, output_filename, verbose=False):

    # Jeroen catalogues compilation
    path_jeroen_compilation = os.path.join(root_dir, '../exoplanets_catalog/table_jeroen')
    name_jeroen_compilation = 'jeroen_cat_compilation.csv'

    # Exoplanet.a catalogue
    path_catalogue_exoplanets_a = os.path.join(root_dir, 'data/exoplanet_data/EXOPLANETS_A')
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
    path_catalogue_exoplanets_org = os.path.join(root_dir, 'data/exoplanet_data/exoplanets.org/')
    name_catalogue_exoplanets_org = 'exoplanets.csv'
    url_catalogue_exoplanets_org = 'http://www.exoplanets.org/csv-files/exoplanets.csv'

    # Nasa Exoplanet Archive catalogue
    path_catalogue_nasa_exoplanet_archive = os.path.join(root_dir, 'data/exoplanet_data/NASAEXOPLANETARCHIVE')
    name_catalogue_nasa_exoplanet_archive = 'nasaexoplanetarchive_full.csv'
    url_catalogue_nasa_exoplanet_archive = ''
    #https://exoplanetarchive.ipac.caltech.edu/workspace/TMP_XTDT5D_5033/TblSearch/2020.06.02_03.25.58_005033/index.html

    # Tepcat allplanet catalogue
    path_allplanet_tepcat = os.path.join(root_dir, 'data/exoplanet_data/tepcat/')
    name_allplanet_tepcat = 'allplanets.csv'
    url_allplanet_tepcat = 'http://www.astro.keele.ac.uk/jkt/tepcat/allplanets-csv.csv'

    # Tepcat observables catalogue
    path_observable_tepcat = os.path.join(root_dir, 'data/exoplanet_data/tepcat/')
    name_observable_tepcat = 'observables.csv'
    url_observable_tepcat = 'http://www.astro.keele.ac.uk/jkt/tepcat/observables.csv'

    # reading catalogues
    jeroen_compilation = get_jeroen_parameters(path_jeroen_compilation, name_jeroen_compilation)

    exoplanet_a_catalogue = get_catalogue(path_catalogue_exoplanets_a,
                                          name_catalogue_exoplanets_a,
                                          url_catalogue=url_catalogue_exoplanets_a, update=False)

    exoplanet_org_catalogue = get_catalogue(path_catalogue_exoplanets_org,
                                            name_catalogue_exoplanets_org,
                                            url_catalogue=url_catalogue_exoplanets_org, update=False)

    nea_catalogue = get_catalogue(path_catalogue_nasa_exoplanet_archive,
                                  name_catalogue_nasa_exoplanet_archive,
                                  url_catalogue=url_catalogue_nasa_exoplanet_archive,
                                  skiprows=303, update=False)

    tepcat_allplanets_catalogue = get_catalogue(path_allplanet_tepcat,
                                                name_allplanet_tepcat,
                                                url_catalogue=url_allplanet_tepcat, update=False)

    tepcat_observable_catalogue = get_catalogue(path_observable_tepcat,
                                                name_observable_tepcat,
                                                url_catalogue=url_observable_tepcat, update=False)

    # getting the list of planets
    list_planets_hst = jeroen_compilation.drop_duplicates('PLANET')['PLANET']

    # reformating catalogues
    exo_a_data = reformat_exoplanets_a_data(exoplanet_a_catalogue)
    exo_org_data = reformat_exoplanets_org_data(exoplanet_org_catalogue)
    nea_data = reformat_nasa_exoplanet_archive_data(nea_catalogue)
    tepcat_data = reformat_tepcat_data(tepcat_allplanets_catalogue, tepcat_observable_catalogue)

    exo_a_merged_hst = merge_catalogue(list_planets_hst, exo_a_data, verbose=False)
    exo_org_merged_hst = merge_catalogue(list_planets_hst, exo_org_data, verbose=False)
    nea_merged_hst = merge_catalogue(list_planets_hst, nea_data, verbose=False)
    tepcat_merged_hst = merge_catalogue(list_planets_hst, tepcat_data, verbose=False)

    # Write merged catalogue with parameters comparison
    #write_diff(exo_a_merged_hst, exo_org_merged_hst, nea_merged_hst, tepcat_merged_hst, output_filename='diff')

    # computing phase error
    exo_a = get_phase_error(jeroen_compilation, exo_a_merged_hst, best_value_only=False)
    exo_org = get_phase_error(jeroen_compilation, exo_org_merged_hst, best_value_only=False)
    nea = get_phase_error(jeroen_compilation, nea_merged_hst, best_value_only=False)
    tepcat = get_phase_error(jeroen_compilation, tepcat_merged_hst, best_value_only=False)

    # Writing results
    writer = pd.ExcelWriter(output_filename + '.xlsx', engine='xlsxwriter')

    jeroen_compilation.to_excel(writer, sheet_name='jeroen_compil', index=False)
    exo_a.to_excel(writer, sheet_name='exo_a', index=False)
    exo_org.to_excel(writer, sheet_name='exo_org', index=False)
    nea.to_excel(writer, sheet_name='nea', index=False)
    tepcat.to_excel(writer, sheet_name='tepcat', index=False)

    exo_a = get_phase_error(jeroen_compilation, exo_a_merged_hst, best_value_only=True)
    exo_org = get_phase_error(jeroen_compilation, exo_org_merged_hst, best_value_only=True)
    nea = get_phase_error(jeroen_compilation, nea_merged_hst, best_value_only=True)
    tepcat = get_phase_error(jeroen_compilation, tepcat_merged_hst, best_value_only=True)
    all_phase_errors = pd.DataFrame({'planet_name': exo_a['planet_name'], 'visit': exo_a['visit'],
                                     'Jeroen': jeroen_compilation['phase_error'],
                                     'exo_a': exo_a['phase_error'], 'exo_org': exo_org['phase_error'],
                                     'nea': nea['phase_error'], 'tepcat': tepcat['phase_error'], 'ttv': nea['ttv']})
    all_phase_errors.to_excel(writer, sheet_name='all_phase_errors', index=False)

    tt_diff = pd.DataFrame({'planet_name': exo_a['planet_name'], 'visit': exo_a['visit'],
                            '(eph_a-eph_j)/eph_j': (exo_a['ephemeris'] - jeroen_compilation['ephemeris']) % jeroen_compilation['period'],
                            '(eph_org-eph_j)/eph_j': (exo_org['ephemeris'] - jeroen_compilation['ephemeris']) % jeroen_compilation['period'],
                            '(eph_nea-eph_j)/eph_j': (nea['ephemeris'] - jeroen_compilation['ephemeris']) % jeroen_compilation['period'],
                            '(eph_tepcat-eph_j)/eph_j': (tepcat['ephemeris'] - jeroen_compilation['ephemeris']) % jeroen_compilation['period'],
                            'ttv': nea['ttv']})
    tt_diff.to_excel(writer, sheet_name='diff_ephemerids', index=False)

    writer.save()


if __name__ == "__main__":
    root_dir = '../'
    # root_dir = os.path.dirname(cascade.__path__[0])   #need to import cascade module
    main(root_dir=root_dir, output_filename='check_cat_phase_error', verbose=False)
