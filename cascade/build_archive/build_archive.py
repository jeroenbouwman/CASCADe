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

@author:Jeroen Bouwman, Rene Gastaud, Raphael Peralta, Fred Lahuis
"""
import os
import stat
import io
import configparser
import shutil
from ast import literal_eval
import numpy as np
import copy
import requests
from tqdm import tqdm
import astropy.units as u
from astropy.table import MaskedColumn
from astropy.io import fits

from ..exoplanet_tools import extract_exoplanet_data
from ..initialize import cascade_default_data_path
from ..initialize import cascade_default_initialization_path
from ..initialize import cascade_default_path
from ..initialize import cascade_default_save_path
from ..initialize import cascade_default_scripts_path
from ..initialize import cascade_default_log_path
from ..initialize import cascade_default_run_scripts_path
from ..initialize import cascade_warnings
from ..initialize import cascade_default_anaconda_environment
from ..exoplanet_tools import parse_database
from ..initialize import read_ini_files


__all__ = ['read_config_file', 'remove_space', 'remove_duplicate',
           'fill_system_parameters', 'long_substr',
           'return_exoplanet_catalogs', 'create_configuration',
           'print_parser_content', 'return_hst_data_calalog_keys',
           'return_all_hst_planets', 'return_header_info',
           'create_bash_script', 'save_observations', 'fill_config_parameters',
           'IniFileParser', 'check_for_exceptions',
           'convert_object_value_strings_to_values', 'create_unique_ids']


def read_config_file(file_name, path):
    """
    Read configuration file to generate .ini file.

    Parameters
    ----------
    file_name : 'str'
        Filename of the configureation file.
    path : 'str'
        Path to the location of the configuration file.

    Returns
    -------
    config_dict : 'dict'
        Dictionary containing the content of the read configuration file.
    """
    config_file = os.path.join(path, file_name)
    with open(config_file, 'r') as conf:
        config_dict = literal_eval(conf.read())
    return config_dict


def remove_space(object_names):
    """
    Remove spaces.

    Parameters
    ----------
    object_names : 'list' or 'ndarray' of 'str'
        List of object names possibly containing spaces.

    Returns
    -------
    new_name : 'ndarray' of 'str'
        List of spaces removed object names.

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
    object_names : 'list' or 'ndarray' of 'str'
        List of object names possibly containing duplicates.

    Returns
    -------
    new_name : 'ndarray' of 'str'
        List of unique object names.

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
        Name of the planetary system.
    catalogs : 'dict'
        Dictionary containing the names of the exoplanet catalogs.
    configuration : 'dict'
        Dictionary containing all relevant entries in the exoplanet catalog.
    primary_cat : 'str'
        primary catalog.

    Raises
    ------
    ValueError
        Errro is raised if system can not be found in catalog.

    Returns
    -------
    system_parameters : 'dict'
        Retrieved exoplanet system parameters from catalog.

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
            data_record = \
                extract_exoplanet_data(catalogs[cat], name,
                                       search_radius=60*u.arcsec)
            no_match = False
        except (ValueError, KeyError) as error:
            print("No match in {}".format(cat))
            print(error)
            no_match = True
            data_record = [{}]
        search_results[cat] = {'no_match': no_match, 'record': data_record}
    # check if system is found in any catalog
    if np.all([i['no_match'] for i in search_results.values()]):
        raise ValueError("Planet not found")
    # does the primary found the target?  if not pick another
    if not search_results[primary_cat]['no_match']:
        data_record = search_results[primary_cat]['record'][0].copy()
    else:
        for cat in search_results:
            if not cat == primary_cat:
                if not search_results[cat]['no_match']:
                    data_record = search_results[cat]['record'][0].copy()
                    break
    # check for missing observable
    for observable, unit in zip(observables, observables_units):
        if observable not in data_record.keys():
            data_record[observable] = MaskedColumn(data=[0.0], mask=True,
                                                   unit=unit)
            for cat in search_results:
                if observable in search_results[cat]['record'][0].keys():
                    data_record[observable] = \
                        search_results[cat]['record'][0][observable]
                    break
    # check for missing values
    for observable in observables:
        if data_record[observable].mask[0]:
            for cat in search_results:
                if observable in search_results[cat]['record'][0].keys():
                    if not search_results[cat]['record'][0][observable].mask[0]:
                        data_record[observable] = \
                             search_results[cat]['record'][0][observable]
                        break
    # Make sure the units are those we use as standard units.
    values = [data_record[key].filled(0.0)[0] *
              data_record[key].unit if key != 'NAME' else
              data_record[key][0] for key in observables]
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
    data : 'ndarray' of 'str'
        Array containing all unique identifiers of the data files.

    Returns
    -------
    substr : 'str'
        Longest common sub-string of the unique identifiers.

    """
    substr = ''
    if len(data) > 1 and len(data[0]) > 0:
        for i in range(len(data[0])):
            for j in range(len(data[0])-i+1):
                if j > len(substr) and all(data[0][i:i+j] in x for x in data):
                    substr = data[0][i:i+j]
    return substr


def return_exoplanet_catalogs(update=True):
    """
    Create dictionary with all exoplanet catalog data.

    Parameters
    ----------
    update : 'bool', optional
    If True update the files csv else use the csv files from local disk.
    Default is True.

    Returns
    -------
    catalog_dict : 'dict'
        Dictionary containing all entries for all exoplanet catalogs.

    """
    all_catalogs = ['TEPCAT', 'EXOPLANETS.ORG', 'NASAEXOPLANETARCHIVE',
                    'EXOPLANETS_A']
    catalog_dict = {}
    for cat in all_catalogs:
        catalog_dict[cat] = parse_database(cat, update=update)
    return catalog_dict


def create_configuration(template, path, parameter_dict):
    """
    Create a parser.

    Parameters
    ----------
    template : 'str'
        DESCRIPTION.
    path : 'str'
        DESCRIPTION.
    parameter_dict : 'dictionary'
        DESCRIPTION.

    Returns
    -------
    parser : 'configparser.ConfigParser'
        DESCRIPTION.

    """
    with open(os.path.join(path, template)) as template_file:
        filled_template = template_file.read().format(**parameter_dict)
        parser = configparser.ConfigParser()
        parser.optionxform = str  # make option names case sensitive
        parser.read_file(io.StringIO(filled_template))
    return parser


def print_parser_content(parser):
    """
    Print parser content.

    Parameters
    ----------
    parser : 'configparser.ConfigParser'
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
    planet_name : 'str'
        DESCRIPTION.
    hst_data_catalog : 'dictionary'
        DESCRIPTION.

    Returns
    -------
    catalog_keys : 'list'
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
    hst_data_catalog : 'dictionary'
        DESCRIPTION.

    Returns
    -------
    all_observed_planets : 'list' of 'str'
        DESCRIPTION.

    """
    all_observed_planets = remove_duplicate(
        [record['planet'] for record in hst_data_catalog.values()]
        )
    return all_observed_planets


def create_unique_ids(hst_data_catalog):
    """
    Check for duplicate ID's.

    Parameters
    ----------
    hst_data_catalog : 'dict'
        DESCRIPTION.

    Returns
    -------
    id_dict_all : 'dict'
        DESCRIPTION.
    """
    numeral = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII']
    id_dict_all = {}
    all_planets = return_all_hst_planets(hst_data_catalog)
    for planet in all_planets:
        visits = return_hst_data_calalog_keys(planet, hst_data_catalog)
        obs_ids = []
        for visit in visits:
            data_files = hst_data_catalog[visit]['observations_id_ima'].split(',')
            data_file_id = [file.split('_')[0] for file in data_files]
            obs_ids.append(long_substr(data_file_id))
        uniq_id = np.unique(obs_ids)
        obs_ids_dir = obs_ids.copy()
        if len(uniq_id) != len(obs_ids):
            for id in uniq_id:
                idx = np.where(np.array(obs_ids) == id)[0]
                if len(idx) > 1:
                    for inum, idx_select in enumerate(idx):
                        obs_ids_dir[idx_select] = \
                            obs_ids_dir[idx_select]+'_'+numeral[inum]
        id_dict_visit = {}
        for visit, obs_id, obs_id_dir in zip(visits, obs_ids, obs_ids_dir):
            id_dict_visit[visit] = {'obs_id': obs_id, 'obs_id_dir': obs_id_dir}
        id_dict_all[planet] = id_dict_visit
    return id_dict_all


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
    fits_keywords = ['TARGNAME', 'RA_TARG', 'DEC_TARG', 'PROPOSID',
                     'TELESCOP', 'INSTRUME', 'FILTER', 'APERTURE',
                     'NSAMP', 'EXPTIME', 'SCAN_TYP', 'SCAN_RAT', 'SCAN_LEN']

    url_data_archive = \
        'https://mast.stsci.edu/portal/Download/file/HST/product/{0}'
    headers = {"User-Agent": "Mozilla/5.0 (X11; CrOS x86_64 12871.102.0) \
               AppleWebKit/537.36 (KHTML, like Gecko) Chrome/81.0.4044.141 \
               Safari/537.36"}

    data_file_id = [file.split('_')[0] for file in [data_file, cal_data_file]]
    data_file_id = long_substr(data_file_id)
    temp_download_dir = os.path.join(cascade_default_data_path,
                                     "mastDownload_"+data_file_id+"/")

    os.makedirs(temp_download_dir, exist_ok=True)
    df = requests.get(url_data_archive.format(data_file), stream=True,
                      headers=headers)
    with open(os.path.join(temp_download_dir, data_file), 'wb') as file:
        for chunk in df.iter_content(chunk_size=1024):
            file.write(chunk)
    header_info = {}
    with fits.open(os.path.join(temp_download_dir, data_file),
                   ignore_missing_end=True) as hdul:
        for keyword in fits_keywords:
            header_info[keyword] = hdul['PRIMARY'].header[keyword]
    df = requests.get(url_data_archive.format(cal_data_file), stream=True,
                      headers=headers)
    with open(os.path.join(temp_download_dir, cal_data_file), 'wb') as file:
        for chunk in df.iter_content(chunk_size=1024):
            file.write(chunk)
    cal_header_info = {}
    with fits.open(os.path.join(temp_download_dir, cal_data_file),
                   ignore_missing_end=True) as hdul:
        for keyword in fits_keywords:
            cal_header_info[keyword] = hdul['PRIMARY'].header[keyword]
    # some house cleaning
    shutil.rmtree(temp_download_dir)

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


def create_bash_script(database_id, configuration):
    """
    Create bash run script.

    Parameters
    ----------
    database_id : 'str'
        DESCRIPTION.
    configuration : 'configparser.ConfigParser'
        DESCRIPTION.

    Returns
    -------
    None.

    """
    # Location of the tamples
    templates_dir = os.path.join(cascade_default_path,
                                 'configuration_templates/')
    scripts_template = 'bash_script.template'

    observatory = configuration['INSTRUMENT']["instrument_observatory"]
    instrument = configuration['INSTRUMENT']["instrument"]
    system_name = configuration['OBSERVATIONS']['observations_target_name']
    instrument_save_path = os.path.join(observatory, instrument, "")
    script_dict = {}
    script_dict['save_path'] = cascade_default_save_path
    script_dict['default_path'] = cascade_default_path
    script_dict['init_path'] = cascade_default_initialization_path
    script_dict['data_path'] = cascade_default_data_path
    script_dict['scripts_path'] = cascade_default_scripts_path 
    script_dict['log_path'] = cascade_default_log_path
    script_dict['warnings'] = cascade_warnings
    script_dict['conda_env'] = cascade_default_anaconda_environment
    script_dict['run_scripts_path'] = cascade_default_run_scripts_path
    script_dict['system_name'] = system_name
    script_dict['instrument_save_path'] = instrument_save_path
    script_dict['visit'] = database_id

    with open(os.path.join(templates_dir, scripts_template)) as template_file:
        filled_template = template_file.read().format(**script_dict)

    #bash_file_path = os.path.join(cascade_default_path, 'scripts',
    #                              instrument_save_path)
    bash_file_path = cascade_default_scripts_path / instrument_save_path
    os.makedirs(bash_file_path, exist_ok=True)
    bash_filename = 'run_'+system_name+'.sh'
    with open(os.path.join(bash_file_path, bash_filename), 'w') as file:
        file.write(filled_template)
    st = os.stat(os.path.join(bash_file_path, bash_filename))
    os.chmod(os.path.join(bash_file_path, bash_filename),
             st.st_mode | stat.S_IEXEC)


def save_observations(data_files, cal_data_files, parser,
                      skip_existing=False):
    """
    Save HST archive data to disk.

    Parameters
    ----------
    data_files : 'list' of 'str'
        DESCRIPTION.
    cal_data_files : 'list' of 'str'
        DESCRIPTION.
    parser : 'configparser.ConfigParser'
        DESCRIPTION.
    skip_excisting : 'bool'
        DESCRIPTION Default is False

    Returns
    -------
    None.

    """
    url_data_archive = \
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

    if parser['OBSERVATIONS']['observations_data'] == 'SPECTRAL_CUBE':
        data_files_to_download = data_files.copy()
    else:
        data_files_to_download = \
            [file.replace('_ima', '_flt') for file in data_files]
    for datafile in tqdm(data_files_to_download, dynamic_ncols=True,
                         desc='Downloading Archive Data '):
        if not (skip_existing and os.path.exists(
                os.path.join(data_save_path, datafile))):
            df = requests.get(url_data_archive.format(datafile), stream=True)
            with open(os.path.join(data_save_path, datafile), 'wb') as file:
                for chunk in df.iter_content(chunk_size=1024):
                    file.write(chunk)
    data_files_to_download = cal_data_files.copy()
    for datafile in tqdm(data_files_to_download, dynamic_ncols=True,
                         desc='Downloading Aquisition Images '):
        if not (skip_existing and os.path.exists(
                os.path.join(data_save_path, datafile))):
            df = requests.get(url_data_archive.format(datafile), stream=True)
            with open(os.path.join(data_save_path, datafile), 'wb') as file:
                for chunk in df.iter_content(chunk_size=1024):
                    file.write(chunk)
    if parser['OBSERVATIONS']['observations_data'] == 'SPECTRAL_CUBE':
        data_files_to_download = \
            [file.replace('_flt', '_ima') for file in cal_data_files]
        for datafile in tqdm(data_files_to_download, dynamic_ncols=True,
                             desc='Downloading Aquisition Images '):
            if not (skip_existing and os.path.exists(
                    os.path.join(data_save_path, datafile))):
                df = requests.get(url_data_archive.format(datafile),
                                  stream=True)
                with open(os.path.join(data_save_path, datafile), 'wb') as file:
                    for chunk in df.iter_content(chunk_size=1024):
                        file.write(chunk)


def fill_config_parameters(config_dict, namespece_dict):
    """
    Fill in values of config dictionary.

    Parameters
    ----------
    config_dict : 'dictionary'
        DESCRIPTION.
    namespece_dict : 'dictionary'
        DESCRIPTION.

    Returns
    -------
    config_dict : 'dictionary'
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
            print("The value {} of the {} parameter does not correspond "
                  "to any of the allowd values: {}. Aboritng creating of "
                  "configuration dictionary".format(new_value, key,
                                                    values['allowed']))
            raise ValueError
    return config_dict


def convert_value_strings_to_values(value_string):
    """
    Convert value strings to values.

    Parameters
    ----------
    value_string : 'str'
        DESCRIPTION.

    Returns
    -------
    value : 'float', 'int', 'bool' or 'str'
        DESCRIPTION.

    """
    try:
        value = literal_eval(value_string)
    except (ValueError, SyntaxError):
        value = value_string
    return value


def convert_object_value_strings_to_values(value_string):
    """
    Convert object value strings to values.

    Parameters
    ----------
    value_string : 'str'
        DESCRIPTION.

    Returns
    -------
    value : 'astropy.units.Quantity' or 'float', 'int', 'bool', 'str'
        DESCRIPTION.

    """
    try:
        value = u.Quantity(value_string)
    except (ValueError, TypeError):
        try:
            value = literal_eval(value_string)
        except (ValueError, SyntaxError):
            value = value_string
    return value


def check_for_exceptions(exception_file, parameter_dict):
    """
    Check for initialization exceptions.

    Parameters
    ----------
    exception_file : 'str'
        DESCRIPTION.
    parameter_dict : 'dictionary'
        DESCRIPTION.

    Returns
    -------
    exceptions_dict : 'dictionary'
        DESCRIPTION.

    """
    observation_name = parameter_dict['observations_target_name']
    observtory = parameter_dict['instrument_observatory']
    instrument = parameter_dict['instrument']
    file_path = os.path.join(cascade_default_path, 'archive_databases/',
                             observtory, instrument)
    if not os.path.isfile(os.path.join(file_path, exception_file)):
        return {}
    parser = read_ini_files(os.path.join(file_path, exception_file))
    section_names = parser.sections()
    if observation_name in section_names:
        exceptions_dict = dict(parser.items(observation_name))
        for key, value_string in exceptions_dict.items():
            if 'object' in key:
                value = convert_object_value_strings_to_values(value_string)
            else:
                value = convert_value_strings_to_values(value_string)
            exceptions_dict[key] = value
    else:
        exceptions_dict = {}
    return exceptions_dict


class IniFileParser:
    """Inititialization file parser class."""

    def __init__(self, configuration_file_list, initialization_file_template,
                 namespace_dict, templates_path, planet_name=None,
                 exoplanet_catalogs_dictionary=None,
                 primary_exoplanet_catalog=None):
        self.configuration_file_list = configuration_file_list
        self.initialization_file_template = initialization_file_template
        self.namespace_dict = namespace_dict
        self.templates_path = templates_path
        self.planet_name = planet_name
        self.exoplanet_catalogs_dictionary = exoplanet_catalogs_dictionary
        self.primary_exoplanet_catalog = primary_exoplanet_catalog
        self.create_parser()

    def create_parser(self):
        """
        Create configuration parser.

        Returns
        -------
        None.

        """
        full_configuration_dict = {}
        for configuration_file in self.configuration_file_list:
            config_dict = \
                read_config_file(configuration_file, self.templates_path)
            if 'object' not in configuration_file:
                config_dict = \
                    fill_config_parameters(config_dict, self.namespace_dict)
            else:
                config_dict = \
                    fill_system_parameters(self.planet_name,
                                           self.exoplanet_catalogs_dictionary,
                                           config_dict,
                                           self.primary_exoplanet_catalog)
                for key in config_dict.keys():
                    if key in self.namespace_dict.keys():
                        config_dict[key] = copy.copy(self.namespace_dict[key])
            full_configuration_dict.update(config_dict.copy())
        init_file_parser = create_configuration(
            self.initialization_file_template,
            self.templates_path,
            full_configuration_dict
            )
        self.init_file_parser = init_file_parser

    def return_parser(self):
        """
        Return configuration parser.

        Returns
        -------
        'configparser.ConfigParser'
            DESCRIPTION.

        """
        return self.init_file_parser

    def print_parser(self):
        """
        Print parser content.

        Returns
        -------
        None.

        """
        print_parser_content(self.init_file_parser)

    def save_parser(self):
        """
        Save parser.

        Returns
        -------
        None.

        """
        pass
