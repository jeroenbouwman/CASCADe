#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# This file is part of the CASCADe package which has been
# developed within the ExoplANETS-A H2020 program.
#
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright (C) 2020  Jeroen Bouwman
"""
Created on April 24 2020

@author:Jeroen Bouwman, Rene Gastaud, Raphael Peralta
"""
import os
import astropy.units as u
import numpy as np
from astropy.table import MaskedColumn
import io
import pickle
import configparser
from urllib.request import urlretrieve
import shutil
from astropy.io import fits
from tqdm import tqdm

os.environ["CASCADE_WARNINGS"] = 'off'
home_dir = os.environ["HOME"]
os.environ['CASCADE_OBSERVATIONS'] = \
    os.path.join(home_dir, 'cascasde_test_hst_data')
os.environ['CASCADE_INITIALIZATION_FILE_PATH'] = \
    os.path.join(home_dir, 'cascade_test_ini_files')

from cascade.initialize import cascade_default_data_path
from cascade.initialize import cascade_default_initialization_path
from cascade.exoplanet_tools import parse_database
from cascade.exoplanet_tools import extract_exoplanet_data


def read_config_file(file_name, path):
    """
    Read configuration file to generate .ini file.

    Parameters
    ----------
    file_name : 'str'
        DESCRIPTION.
    path : 'str'
        DESCRIPTION.

    Returns
    -------
    config_dict : 'dict'
        DESCRIPTION.
    """
    from ast import literal_eval
    config_file = os.path.join(path, file_name)
    with open(config_file, 'r') as conf:
        config_dict = literal_eval(conf.read())
    return config_dict


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


def fill_system_parameters(name, catalogs, configuration,
                           primary_cat):
    """
    Return observed system parameters.

    Parameters
    ----------
    name : 'str'
        DESCRIPTION.
    catalogs : 'dict'
        DESCRIPTION.
    configuration : 'dict'
        DESCRIPTION
    primary_cat : 'str'
        DESCRIPTION.

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    system_parameters : 'dict'
        DESCRIPTION.

    """
    observables = []
    observables_units = []
    parameters = []
    for parameter, catalog_entry in configuration.items():
        parameters.append(parameter)
        observables.append(catalog_entry['key'])
        observables_units.append(u.Unit(catalog_entry['unit']))

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
    values = [dr[key].filled(0.0)[0] * dr[key].unit if key != 'NAME' else
              dr[key][0] for key in observables]
    for i, (value, unit) in enumerate(zip(values, observables_units)):
        if unit != u.dimensionless_unscaled:
            values[i] = value.to(unit)
    system_parameters = dict(zip(parameters, values))
    return system_parameters


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


def create_configuration(template, path, parameter_dict):
    """
    Create a parser.

    Parameters
    ----------
    template : TYPE
        DESCRIPTION.
    path : TYPE
        DESCRIPTION.
    parameter_dict : TYPE
        DESCRIPTION.

    Returns
    -------
    parser : TYPE
        DESCRIPTION.

    """
    with open(os.path.join(path, template)) as template_file:
        filled_template = template_file.read().format(**parameter_dict)
        parser = configparser.ConfigParser()
        parser.optionxform = str  # make option names case sensitive
        parser.readfp(io.StringIO(filled_template))
    return parser


def print_parser_content(parser):
    """
    Print parser content.

    Parameters
    ----------
    parser : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    for section_name in parser.sections():
        print('[{}]'.format(section_name))
        for parameter, value in parser.items(section_name):
            print('%s = %s' % (parameter, value))
        print()


def return_hst_data_calalog_keys(planet_name, hst_data_catalog):
    """
    Return catalog keys for planet.

    Parameters
    ----------
    planet_name : TYPE
        DESCRIPTION.
    hst_data_catalog : TYPE
        DESCRIPTION.

    Returns
    -------
    catalog_keys : TYPE
        DESCRIPTION.

    """
    catalog_keys = [key for key, record in hst_data_catalog.items()
                    if record['planet'] == planet_name]

    return catalog_keys


def return_all_hst_planets(hst_data_catalog):
    """
    Return list with all observed planets.

    Parameters
    ----------
    hst_data_catalog : TYPE
        DESCRIPTION.

    Returns
    -------
    all_observed_planets : TYPE
        DESCRIPTION.

    """
    all_observed_planets = remove_duplicate(
        [record['planet'] for record in hst_data_catalog.values()]
        )
    return all_observed_planets


def return_header_info(data_file, cal_data_file):
    """
    Return all relavant parameters from fits data file header.

    Parameters
    ----------
    data_file : 'str'
        DESCRIPTION.
    cal_data_file : 'str'
        DESCRIPTION.

    Returns
    -------
    cascade_parameter_dict : 'dict'
        DESCRIPTION.

    """
    FITS_KEYWORDS = ['TARGNAME', 'RA_TARG', 'DEC_TARG', 'PROPOSID',
                     'TELESCOP', 'INSTRUME', 'FILTER', 'APERTURE',
                     'NSAMP', 'EXPTIME', 'SCAN_TYP', 'SCAN_RAT', 'SCAN_LEN']

    URL_DATA_ARCHIVE = \
        'https://mast.stsci.edu/portal/Download/file/HST/product/{0}'

    TEMP_DOWNLOAD_DIR = os.path.join(cascade_default_data_path,
                                     "mastDownload/")

    os.makedirs(TEMP_DOWNLOAD_DIR, exist_ok=True)
    urlretrieve(URL_DATA_ARCHIVE.format(data_file),
                os.path.join(TEMP_DOWNLOAD_DIR, data_file))
    header_info = {}
    with fits.open(os.path.join(TEMP_DOWNLOAD_DIR, data_file)) as hdul:
        for keyword in FITS_KEYWORDS:
            header_info[keyword] = hdul['PRIMARY'].header[keyword]

    urlretrieve(URL_DATA_ARCHIVE.format(cal_data_file),
                os.path.join(TEMP_DOWNLOAD_DIR, cal_data_file))
    cal_header_info = {}
    with fits.open(os.path.join(TEMP_DOWNLOAD_DIR, cal_data_file)) as hdul:
        for keyword in FITS_KEYWORDS:
            cal_header_info[keyword] = hdul['PRIMARY'].header[keyword]
    # some house cleaning
    shutil.rmtree(TEMP_DOWNLOAD_DIR)

    if header_info['SCAN_TYP'] == 'N':
        obs_mode = 'STARING'
        obs_data = 'SPECTRAL_IMAGE'
        data_product = 'flt'
        nextraction = 7
    else:
        obs_mode = 'SCANNING'
        obs_data = 'SPECTRAL_CUBE'
        data_product = 'ima'
        pixel_sze = 0.121
        number_of_samples = int(header_info['NSAMP'])-2
        scan_length = float(header_info['SCAN_LEN'])
        nextraction = int(scan_length/pixel_sze) // number_of_samples + 7
        if nextraction % 2 == 0:  # even
            nextraction += 1
        nextraction = str(nextraction)

    cascade_parameter_dict = {}
    cascade_parameter_dict['processing_nextraction'] = nextraction
    cascade_parameter_dict['observations_mode'] = obs_mode
    cascade_parameter_dict['observations_data'] = obs_data
    cascade_parameter_dict['observations_data_product'] = data_product
    cascade_parameter_dict['observations_has_background'] = True
    cascade_parameter_dict['instrument_observatory'] = header_info['TELESCOP']
    cascade_parameter_dict['instrument'] = header_info['INSTRUME']
    cascade_parameter_dict['instrument_filter'] = header_info['FILTER']
    cascade_parameter_dict['instrument_aperture'] = header_info['APERTURE']
    cascade_parameter_dict['instrument_cal_filter'] = cal_header_info['FILTER']
    cascade_parameter_dict['instrument_cal_aperture'] = \
        cal_header_info['APERTURE']
    return cascade_parameter_dict


def save_observations(data_files, cal_data_files, parser):
    """
    Save HST archive data to disk.

    Parameters
    ----------
    data_files : TYPE
        DESCRIPTION.
    cal_data_files : TYPE
        DESCRIPTION.
    parser : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    URL_DATA_ARCHIVE = \
        'https://mast.stsci.edu/portal/Download/file/HST/product/{0}'

    data_save_path = os.path.join(
        cascade_default_data_path,
        parser['OBSERVATIONS']['observations_path'],
        parser['INSTRUMENT']["instrument_observatory"],
        parser['INSTRUMENT']["instrument"],
        parser['OBSERVATIONS']['observations_target_name'],
        "SPECTRAL_IMAGES"
        )
    os.makedirs(data_save_path, exist_ok=True)

    if parser['OBSERVATIONS']['observations_mode'] == 'SCANNING':
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
    if parser['OBSERVATIONS']['observations_mode'] == 'SCANNING':
        data_files_to_download = \
            [file.replace('_flt', '_ima') for file in cal_data_files]
        for datafile in tqdm(data_files_to_download, dynamic_ncols=True):
            urlretrieve(URL_DATA_ARCHIVE.format(datafile),
                        os.path.join(data_save_path, datafile))


def fill_config_parameters(config_dict, namespece_dict):
    """
    Fill in values of config dictionary.

    Parameters
    ----------
    config_dict : TYPE
        DESCRIPTION.
    namespece_dict : TYPE
        DESCRIPTION.

    Returns
    -------
    config_dict : TYPE
        DESCRIPTION.

    """
    for key, values in config_dict.items():
        new_value = namespece_dict.get(key, values['default'])
        if new_value == 'NO_CASCADE_DEFAULT_VALUE':
            print("The {} parameter has no defaut value and is not "
                  "defined. Aborting creaton of configuration "
                  "dictionary".format(key))
            raise ValueError
        if ((new_value in values['allowed']) |
                (values['allowed'][0] == 'NO_RESTRICTIONS')):
            config_dict[key] = new_value
        else:
            print("The value of the {} parameter does not correspond "
                  "to any of the allowd values: {}. Aboritng creating of "
                  "configuration dictionary".format(key, values['allowed']))
            raise ValueError
    return config_dict


# ####################### START ###########################

# Primary exoplanet catalog
PRIMARY_CATALOG = "NASAEXOPLANETARCHIVE"

# Location of the tamples
TEMPLATES_DIR = os.path.join(cascade_default_data_path,
                             'configuration_templates/')

# All configuration files for the different sections in the
# initialization files
OBJECT_CONFIGURATION_FILE = 'object.conf'
CATALOG_CONFIGURATION_FILE = 'catalog.conf'
CASCADE_CONFIGURATION_FILE = 'cascade.conf'
PROCESSING_CONFIGURATION_FILE = 'processing.conf'
INSTRUMENT_CONFIGURATION_FILE = 'instrument.conf'
OBSERVATIONS_CONFIGURATION_FILE = 'observations.conf'
CPM_CONFIGURATION_FILE = 'cpm.conf'
MODEL_CONFIGURATION_FILE = 'model.conf'

# All initialization filetemplates
OBJECT_CONFIGURATION_TEMPLATE = 'cascade_object_template.ini'
EXTRACT_TIMESERIES_CONFIGURATION_TEMPLATE = \
    'cascade_extract_timeseries_template.ini'
CALIBRATE_PLANET_SPECTRUM_CONFIGURATION_TEMPLATE = \
    'cascade_calibrate_planet_spectrum_template.ini'

# Get the HST observations catalog file
HST_CATALOG_FILE = os.path.join(cascade_default_data_path,
                                "archive_databases/HST/WFC3/",
                                "WFC3_files.pickle")
with open(HST_CATALOG_FILE, 'rb') as f:
    hst_data_catalog = pickle.load(f)

# explanet data catalogs
catalogs_dictionary = return_exoplanet_catalogs()


# TEST
print(return_all_hst_planets(hst_data_catalog))
print(return_hst_data_calalog_keys('WASP-19b', hst_data_catalog))
visits = ['12181.14']

# ########### HERE LOOP STARTS #####################
if (visits is None):
    visits = hst_data_catalog.keys()
for visit in visits:
    PLANET_NAME = hst_data_catalog[visit]['planet']
    print("Target: {}, Visit: {}".format(PLANET_NAME, visit))

    # ########3####  some logic  ##############
    # ## Define the observation typpe and skips loop if problem ##
    if hst_data_catalog[visit]['observation'] == 'transit':
        OBS_TYPE = 'TRANSIT'
    elif hst_data_catalog[visit]['observation'] == 'eclipse':
        OBS_TYPE = 'ECLIPSE'
    else:
        OBS_TYPE = 'PROBLEM'
    if OBS_TYPE == 'PROBLEM':
        print('Skipping : {}'.format(hst_data_catalog[visit]['observation']))
        continue

    # ################ CREATE OBJECT.INI ##############################
    # read object.conf file
    object_config_dict = \
        read_config_file(OBJECT_CONFIGURATION_FILE, TEMPLATES_DIR)
    # fill with exoplanet catalog data
    object_config_dict = \
        fill_system_parameters(PLANET_NAME, catalogs_dictionary,
                               object_config_dict, PRIMARY_CATALOG)
    # read catalog.conf file and fill parameters with values.
    # set which catalog to use
    catalog_name = PRIMARY_CATALOG
    # pass along namespace dict containing relevant parameters
    namespace_dict = {**locals()}
    catalog_config_dict = \
        read_config_file(CATALOG_CONFIGURATION_FILE, TEMPLATES_DIR)
    catalog_config_dict = \
        fill_config_parameters(catalog_config_dict, namespace_dict)

    # combine configureation
    combined_object_initialization_file_dict = \
        {**object_config_dict,
         **catalog_config_dict
         }
    # create parser
    object_parser = create_configuration(
        OBJECT_CONFIGURATION_TEMPLATE,
        TEMPLATES_DIR,
        combined_object_initialization_file_dict)
    print_parser_content(object_parser)

    # ################ CREATE EXTRACT_TIMESERIES.INI ####################
    create_timeseries_namespace_dict = {}
    # get the data files from hst databes
    data_files = hst_data_catalog[visit]['observations_id_ima'].split(',')
    cal_data_files = hst_data_catalog[visit]['calibrations_id'].split(',')
    data_file_id = [file.split('_')[0] for file in data_files]
    cal_data_file_id = [file.split('_')[0] for file in cal_data_files]

    # get or calculate all parameters needed for the observations and
    # instrument sections in the .ini file
    OBS_ID = long_substr(data_file_id)

    # fill dictionary with parameters to be filled into tamplates
    create_timeseries_namespace_dict['cascade_save_path'] = \
        PLANET_NAME+'_'+OBS_ID
    create_timeseries_namespace_dict['observations_type'] = OBS_TYPE
    create_timeseries_namespace_dict['observations_target_name'] = \
        PLANET_NAME+'_'+OBS_ID
    create_timeseries_namespace_dict['observations_id'] = OBS_ID
    create_timeseries_namespace_dict.update(
        **return_header_info(data_files[0], cal_data_files[0])
        )

    # get the .conf files and fill with values.
    cascade_config_dict = \
        read_config_file(CASCADE_CONFIGURATION_FILE, TEMPLATES_DIR)
    cascade_config_dict = \
        fill_config_parameters(cascade_config_dict,
                               create_timeseries_namespace_dict)
    processing_config_dict = \
        read_config_file(PROCESSING_CONFIGURATION_FILE, TEMPLATES_DIR)
    processing_config_dict = \
        fill_config_parameters(processing_config_dict,
                               create_timeseries_namespace_dict)
    instrument_config_dict = \
        read_config_file(INSTRUMENT_CONFIGURATION_FILE, TEMPLATES_DIR)
    instrument_config_dict = \
        fill_config_parameters(instrument_config_dict,
                               create_timeseries_namespace_dict)
    observations_config_dict = \
        read_config_file(OBSERVATIONS_CONFIGURATION_FILE, TEMPLATES_DIR)
    observations_config_dict = \
        fill_config_parameters(observations_config_dict,
                               create_timeseries_namespace_dict)

    # combine configuration
    combined_extract_timeseries_initialization_file_dict = \
        {**cascade_config_dict,
         **processing_config_dict,
         **instrument_config_dict,
         **observations_config_dict
         }
    # create parser
    extract_timeseries_parser = create_configuration(
        EXTRACT_TIMESERIES_CONFIGURATION_TEMPLATE,
        TEMPLATES_DIR,
        combined_extract_timeseries_initialization_file_dict
        )
    print_parser_content(extract_timeseries_parser)

    # ################ CREATE CALIBRATE_PLANET_SPECTRUM.INI ##################
    # update parameters for timeseries of 1D spectra
    calibrate_planet_spectrum_namespace_dict = \
        create_timeseries_namespace_dict.copy()
    calibrate_planet_spectrum_namespace_dict['processing_nextraction'] = 1
    calibrate_planet_spectrum_namespace_dict['observations_data_product'] = \
        'COE'
    calibrate_planet_spectrum_namespace_dict['observations_has_background'] = \
        False
    calibrate_planet_spectrum_namespace_dict['observations_mode'] = 'STARING'
    calibrate_planet_spectrum_namespace_dict['observations_data'] = 'SPECTRUM'

    # read the .conf and fill values
    processing_config_dict = \
        read_config_file(PROCESSING_CONFIGURATION_FILE, TEMPLATES_DIR)
    processing_config_dict = \
        fill_config_parameters(processing_config_dict,
                               calibrate_planet_spectrum_namespace_dict)
    cpm_config_dict = \
        read_config_file(CPM_CONFIGURATION_FILE, TEMPLATES_DIR)
    cpm_config_dict = \
        fill_config_parameters(cpm_config_dict,
                               calibrate_planet_spectrum_namespace_dict)
    model_config_dict = \
        read_config_file(MODEL_CONFIGURATION_FILE, TEMPLATES_DIR)
    model_config_dict = \
        fill_config_parameters(model_config_dict,
                               calibrate_planet_spectrum_namespace_dict)
    observations_config_dict = \
        read_config_file(OBSERVATIONS_CONFIGURATION_FILE, TEMPLATES_DIR)
    observations_config_dict = \
        fill_config_parameters(observations_config_dict,
                               calibrate_planet_spectrum_namespace_dict)
    # combine configuration
    combined_calibrate_planet_spectrum_initialization_file_dict = \
        {**cascade_config_dict,
         **processing_config_dict,
         **cpm_config_dict,
         **model_config_dict,
         **instrument_config_dict,
         **observations_config_dict
         }
    # create parser
    calibrate_planet_spectrum_parser = create_configuration(
            CALIBRATE_PLANET_SPECTRUM_CONFIGURATION_TEMPLATE,
            TEMPLATES_DIR,
            combined_calibrate_planet_spectrum_initialization_file_dict
            )
    print_parser_content(calibrate_planet_spectrum_parser)

    # ############## Saving ini files ##################
    initialization_save_path = os.path.join(
        cascade_default_initialization_path,
        extract_timeseries_parser['INSTRUMENT']["instrument_observatory"],
        extract_timeseries_parser['INSTRUMENT']["instrument"],
        extract_timeseries_parser['OBSERVATIONS']['observations_target_name']
        )
    os.makedirs(initialization_save_path, exist_ok=True)

    configuration_base_filename = \
        (
         "cascade_" +
         extract_timeseries_parser['OBSERVATIONS']
         ['observations_target_name']
         )
    configuration_object_filename = \
        configuration_base_filename+"_object.ini"
    with open(os.path.join(initialization_save_path,
                           configuration_object_filename), 'w') as configfile:
        object_parser.write(configfile)

    configuration_extract_timeseries_filename = \
        configuration_base_filename+"_extract_timeseries.ini"
    with open(os.path.join(
            initialization_save_path,
            configuration_extract_timeseries_filename), 'w') as configfile:
        extract_timeseries_parser.write(configfile)

    configuration_calibrate_planet_spectrum_filename = \
        configuration_base_filename+"_calibrate_planet_spectrum.ini"
    with open(os.path.join(
            initialization_save_path,
            configuration_calibrate_planet_spectrum_filename),
            'w') as configfile:
        calibrate_planet_spectrum_parser.write(configfile)

    # # ################# GET ARCHIVE DATA ######################
    save_observations(data_files, cal_data_files, extract_timeseries_parser)
