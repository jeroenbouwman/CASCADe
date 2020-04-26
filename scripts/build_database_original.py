#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 09:16:41 2020

@author: bouwman
"""
from urllib.request import urlretrieve
import pickle
import os
import astropy.units as u
import numpy as np
from astropy.table import MaskedColumn
import shutil
from astropy.io import fits
import configparser
from tqdm import tqdm

os.environ["CASCADE_WARNINGS"] = 'off'

home_dir = os.environ["HOME"]
os.environ['CASCADE_DATA_PATH'] = \
    os.path.join(home_dir, 'cascasde_test_hst_data')
os.environ['CASCADE_INITIALIZATION_FILE_PATH'] = \
    os.path.join(home_dir, 'cascade_test_ini_files')

from cascade.exoplanet_tools import Kmag
from cascade import __path__
from cascade.initialize import cascade_default_data_path
from cascade.initialize import cascade_default_initialization_path
from cascade.exoplanet_tools import parse_database
from cascade.exoplanet_tools import extract_exoplanet_data


def remove_space(object_names):
    """
    Remove spaces.

    Parameters
    ----------
    object_names : TYPE
        DESCRIPTION.

    Returns
    -------
    new_name : TYPE
        DESCRIPTION.

    """
    new_name = np.asarray(object_names)
    for i in range(new_name.size):
        new_name[i] = (''.join(object_names[i].split(' ')))
    return new_name


def remove_duplicate(object_names):
    """
    Remove duplicates.

    Parameters
    ----------
    object_names : TYPE
        DESCRIPTION.

    Returns
    -------
    new_name : TYPE
        DESCRIPTION.

    """
    new_name = np.asarray(object_names)
    new_name = np.unique(new_name)
    return new_name


def return_system_parameters(name, catalogs, observables, observables_units,
                             primary_cat):
    """
    Return observed system parameters.

    Parameters
    ----------
    name : TYPE
        DESCRIPTION.
    catalogs : TYPE
        DESCRIPTION.
    observables : TYPE
        DESCRIPTION.
    observables_units : TYPE
        DESCRIPTION.
    primary_cat : TYPE
        DESCRIPTION.

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    values : TYPE
        DESCRIPTION.

    """
    search_results = {}
    for cat in catalogs:
        try:
            dr = extract_exoplanet_data(catalogs[cat],
                                        name, search_radius=2*u.arcsec)
            no_match = False
        except (ValueError, KeyError) as e:
            print("No match in {}".format(cat))
            print(e)
            no_match = True
            dr = [{}]
        search_results[cat] = {'no_match': no_match, 'record': dr}
    # check if system is found in any catalog
    if np.all([i['no_match'] for i in search_results.values()]):
        raise ValueError("Planet not found")
    # does the primary found the target?  if not pick another
    if not search_results[primary_cat]['no_match']:
        dr = search_results[primary_cat]['record'][0].copy()
    else:
        for cat in search_results:
            if not cat == primary_cat:
                if not search_results[cat]['no_match']:
                    dr = search_results[cat]['record'][0].copy()
                    break
    # check for missing observable
    for observable, unit in zip(observables, observables_units):
        if observable not in dr.keys():
            dr[observable] = MaskedColumn(data=[0.0], mask=True, unit=unit)
            for cat in search_results:
                if observable in search_results[cat]['record'][0].keys():
                    dr[observable] = \
                        search_results[cat]['record'][0][observable]
                    break
    # check for missing values
    for observable in observables:
        if dr[observable].mask[0]:
            for cat in search_results:
                if observable in search_results[cat]['record'][0].keys():
                    if not search_results[cat]['record'][0][observable].mask[0]:
                        dr[observable] = \
                             search_results[cat]['record'][0][observable]
                        break
    # Make sure the units are those we use as standard units.
    values = [dr[key].filled(0.0)[0] * dr[key].unit for key in observables]
    for i, (value, unit) in enumerate(zip(values, observables_units)):
        values[i] = value.to(unit)
    values = tuple(values)
    return values


def long_substr(data):
    """
    Find longest common substring.

    Parameters
    ----------
    data : TYPE
        DESCRIPTION.

    Returns
    -------
    substr : TYPE
        DESCRIPTION.

    """
    substr = ''
    if len(data) > 1 and len(data[0]) > 0:
        for i in range(len(data[0])):
            for j in range(len(data[0])-i+1):
                if j > len(substr) and all(data[0][i:i+j] in x for x in data):
                    substr = data[0][i:i+j]
    return substr


def return_exoplanet_catalogs():
    """
    Create dictionary with all exoplanet catalog data.

    Returns
    -------
    catalog_dict : 'dict'
        DESCRIPTION.

    """
    all_catalogs = ['TEPCAT', 'EXOPLANETS.ORG', 'NASAEXOPLANETARCHIVE',
                    'EXOPLANETS_A']
    catalog_dict = {}
    for cat in all_catalogs:
        catalog_dict[cat] = parse_database(cat, update=True)
    return catalog_dict


def add_configuration(config, section_name, parameter_dict):
    """
    Add definitions to the configuration file.

    Parameters
    ----------
    config : 'configparser'
        DESCRIPTION.
    section_name : 'str'
        DESCRIPTION
    obs_dict : 'dict'
        DESCRIPTION.

    Returns
    -------
    None
    """
    configurations = {}
    for key, value in parameter_dict.items():
        configurations[key] = value
        config[section_name] = configurations


PRIMARY_CATALOG = "NASAEXOPLANETARCHIVE"

OBSERVABLES = ['TT', 'PER', 'R', 'A', 'I', 'ECC', 'OM',
               'RSTAR', 'TEFF', 'KMAG', 'LOGG', 'FE']

OBSERVABLES_UNITS = [u.day, u.day, u.Rjup, u.AU, u.deg,
                     u.dimensionless_unscaled, u.deg,
                     u.Rsun, u.K, Kmag, u.dex(u.cm/u.s**2), u.dex]

FITS_KEYWORDS = ['TARGNAME', 'RA_TARG', 'DEC_TARG', 'PROPOSID',
                 'TELESCOP', 'INSTRUME', 'FILTER', 'APERTURE',
                 'NSAMP', 'EXPTIME', 'SCAN_TYP', 'SCAN_RAT', 'SCAN_LEN']

URL_DATA_ARCHIVE = \
    'https://mast.stsci.edu/portal/Download/file/HST/product/{0}'

HST_CATALOG_FILE = os.path.join(os.path.dirname(__path__[0]),
                                "data/archive_databases/HST/WFC3/",
                                "WFC3_files.pickle")

cat_dict = return_exoplanet_catalogs()

with open(HST_CATALOG_FILE, 'rb') as f:
    hst_data = pickle.load(f)

all_observed_planets = [record['planet'] for record in hst_data.values()]
all_observed_planets = remove_duplicate(all_observed_planets)


temp_dir_path = os.path.join(cascade_default_data_path,
                             "mastDownload/")


observations_definition_dict = {}
instrument_definition_dict = {}
cascade_definitions_dict = {}
processing_definitions_dict = {}
object_definition_dict = {}
catalog_definition_dict = {}
config = configparser.ConfigParser()
config_object = configparser.ConfigParser()
config.optionxform = str
config_object.optionxform = str

#hst_data = {'14260.16': hst_data['14260.16'].copy()}
hst_data = {'15135.02': hst_data['15135.02'].copy()}
# loop over all entries in HST observations catalog
for visit in hst_data:
    name = hst_data[visit]['planet']
    print("Target: {}, Visit: {}".format(name, visit))
    try:
        system_parameters = \
            return_system_parameters(name, cat_dict, OBSERVABLES,
                                     OBSERVABLES_UNITS, PRIMARY_CATALOG)
    except ValueError as e:
        print(e)
        continue
    object_dict = {}
    for key, value in zip(OBSERVABLES, system_parameters):
        object_dict[key] = value

    object_definition_dict['object_name'] = name
    object_definition_dict['object_ephemeris'] = object_dict['TT']
    object_definition_dict['object_period'] = object_dict['PER']
    object_definition_dict['object_radius'] = object_dict['R']
    object_definition_dict['object_semi_major_axis'] = object_dict['A']
    object_definition_dict['object_inclination'] = object_dict['I']
    object_definition_dict['object_eccentricity'] = object_dict['ECC']
    object_definition_dict['object_omega'] = object_dict['OM']
    object_definition_dict['object_radius_host_star'] = object_dict['RSTAR']
    object_definition_dict['object_temperature_host_star'] = \
        object_dict['TEFF']
    object_definition_dict['object_kmag'] = object_dict['KMAG']
    object_definition_dict['object_metallicity_host_star'] = object_dict['FE']
    object_definition_dict['object_logg_host_star'] = object_dict['LOGG']
    add_configuration(config_object, 'OBJECT', object_definition_dict)

    catalog_definition_dict['catalog_use_catalog'] = 'False'
    catalog_definition_dict['catalog_name'] = PRIMARY_CATALOG
    catalog_definition_dict['catalog_update'] = 'True'
    add_configuration(config_object, 'CATALOG', catalog_definition_dict)

    if hst_data[visit]['observation'] == 'transit':
        observations_type = 'TRANSIT'
    elif hst_data[visit]['observation'] == 'eclipse':
        observations_type = 'ECLIPSE'
    else:
        observations_type = 'PROBLEM'

    if observations_type == 'PROBLEM':
        print('Scipping problem: {}'.format(hst_data[visit]['observation']))
        continue
    else:
        observations_definition_dict['observations_type'] = observations_type

    data_files = hst_data[visit]['observations_id_ima'].split(',')
    cal_data_files = hst_data[visit]['calibrations_id'].split(',')
    obsid = [file.split('_')[0] for file in data_files]
    cal_obsid = [file.split('_')[0] for file in cal_data_files]

    common_id = long_substr(obsid)
    # print("Common ID: {} \n".format(common_id))

    # check spectral images
    os.makedirs(temp_dir_path, exist_ok=True)
    urlretrieve(URL_DATA_ARCHIVE.format(data_files[0]),
                os.path.join(temp_dir_path, data_files[0]))
    header_info = {}
    with fits.open(os.path.join(temp_dir_path, data_files[0])) as hdul:
        for keyword in FITS_KEYWORDS:
            header_info[keyword] = hdul['PRIMARY'].header[keyword]

    if header_info['SCAN_TYP'] == 'N':
        observations_mode = 'STARING'
        observations_data = 'SPECTRAL_IMAGE'
        observations_data_product = 'flt'
        processing_nextraction = '7'
    else:
        observations_mode = 'SCANNING'
        observations_data = 'SPECTRAL_CUBE'
        observations_data_product = 'ima'
        pixel_sze = 0.121
        number_of_samples = int(header_info['NSAMP'])-2
        scan_length = float(header_info['SCAN_LEN'])
        nextraction = int(scan_length/pixel_sze) // number_of_samples + 7
        if nextraction % 2 == 0:  # even
            nextraction += 1
        processing_nextraction = str(nextraction)

    cascade_definitions_dict['cascade_save_path'] = name+'_'+common_id
    cascade_definitions_dict['cascade_use_multi_processes'] = 'True'
    cascade_definitions_dict['cascade_verbose'] = 'True'
    cascade_definitions_dict['cascade_save_verbose'] = 'True'
    add_configuration(config, 'CASCADE', cascade_definitions_dict)

    processing_definitions_dict['processing_sigma_filtering'] = '3.5'
    processing_definitions_dict['processing_max_number_of_iterations_filtering'] = '15'
    processing_definitions_dict['processing_fractional_acceptance_limit_filtering'] = '0.005'
    processing_definitions_dict['processing_quantile_cut_movement'] = '0.1'
    processing_definitions_dict['processing_order_trace_movement'] = '1'
    processing_definitions_dict['processing_nreferences_movement'] = '6'
    processing_definitions_dict['processing_main_reference_movement'] = '4'
    processing_definitions_dict['processing_upsample_factor_movement'] = '111'
    processing_definitions_dict['processing_angle_oversampling_movement'] = '2'
    processing_definitions_dict['processing_nextraction'] = \
        processing_nextraction
    processing_definitions_dict['processing_rebin_factor_extract1d'] = '1.05'
    add_configuration(config, 'PROCESSING', processing_definitions_dict)

    observations_definition_dict['observations_mode'] = observations_mode
    observations_definition_dict['observations_data'] = observations_data
    observations_definition_dict['observations_path'] = 'data/HST/'
    observations_definition_dict['observations_target_name'] = \
        name+'_'+common_id
    observations_definition_dict['observations_cal_path'] = 'calibration/HST/'
    observations_definition_dict['observations_id'] = common_id
    observations_definition_dict['observations_cal_version'] = '4.32'
    observations_definition_dict['observations_data_product'] = \
        observations_data_product
    observations_definition_dict['observations_has_background'] = 'True'
    observations_definition_dict['observations_uses_background_model'] = 'True'
    observations_definition_dict['observations_median_signal'] = '0.01'
    add_configuration(config, 'OBSERVATIONS', observations_definition_dict)

    # check aquisition images
    urlretrieve(URL_DATA_ARCHIVE.format(cal_data_files[0]),
                os.path.join(temp_dir_path, cal_data_files[0]))
    cal_header_info = {}
    with fits.open(os.path.join(temp_dir_path, cal_data_files[0])) as hdul:
        for keyword in FITS_KEYWORDS:
            cal_header_info[keyword] = hdul['PRIMARY'].header[keyword]

    instrument_definition_dict['instrument_observatory'] = \
        header_info['TELESCOP']
    instrument_definition_dict['instrument'] = header_info['INSTRUME']
    instrument_definition_dict['instrument_filter'] = header_info['FILTER']
    instrument_definition_dict['instrument_aperture'] = header_info['APERTURE']
    instrument_definition_dict['instrument_cal_filter'] = \
        cal_header_info['FILTER']
    instrument_definition_dict['instrument_cal_aperture'] = \
        cal_header_info['APERTURE']
    instrument_definition_dict['instrument_beam'] = 'A'
    add_configuration(config, 'INSTRUMENTS', instrument_definition_dict)

    shutil.rmtree(os.path.join(temp_dir_path))

    # print and save configuration files
    for section_name in config_object.sections():
        print('[{}]'.format(section_name))
        for parameter, value in config_object.items(section_name):
            print('%s = %s' % (parameter, value))
        print()
    for section_name in config.sections():
        print('[{}]'.format(section_name))
        for parameter, value in config.items(section_name):
            print('%s = %s' % (parameter, value))
        print()
    configuration_save_path = \
        os.path.join(cascade_default_initialization_path,
                     instrument_definition_dict['instrument_observatory'],
                     instrument_definition_dict['instrument'],
                     name+'_'+common_id)
    os.makedirs(configuration_save_path, exist_ok=True)
    configuration_datafile_name = \
        "cascade_"+name+'_'+common_id+"_extract_spectra.ini"
    with open(os.path.join(configuration_save_path,
                           configuration_datafile_name), 'w') as configfile:
        config.write(configfile)
    configuration_objectfile_name = \
        "cascade_"+name+'_'+common_id+"_object.ini"
    with open(os.path.join(configuration_save_path,
                           configuration_objectfile_name), 'w') as configfile:
        config_object.write(configfile)

    # download all data
    data_save_path = \
        os.path.join(cascade_default_data_path,
                     instrument_definition_dict['instrument_observatory'],
                     instrument_definition_dict['instrument'],
                     name+'_'+common_id, "SPECTRAL_IMAGES")
    os.makedirs(data_save_path, exist_ok=True)
    if observations_definition_dict['observations_mode'] == 'SCANNING':
        data_files_to_download = data_files.copy()
    else:
        data_files_to_download = \
            [file.replace('_ima', '_flt') for file in data_files]
    for datafile in tqdm(data_files_to_download, dynamic_ncols=True):
        urlretrieve(URL_DATA_ARCHIVE.format(datafile),
                    os.path.join(data_save_path, datafile))
    data_files_to_download = cal_data_files.copy()
    for datafile in tqdm(data_files_to_download, dynamic_ncols=True):
        urlretrieve(URL_DATA_ARCHIVE.format(datafile),
                    os.path.join(data_save_path, datafile))
    if observations_definition_dict['observations_mode'] == 'SCANNING':
        data_files_to_download = \
            [file.replace('_flt', '_ima') for file in cal_data_files]
        for datafile in tqdm(data_files_to_download, dynamic_ncols=True):
            urlretrieve(URL_DATA_ARCHIVE.format(datafile),
                        os.path.join(data_save_path, datafile))
    break

