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
import click
import sys
import six
from pyfiglet import figlet_format

# os.environ["CASCADE_WARNINGS"] = 'off'
# home_dir = os.environ["HOME"]
# os.environ['CASCADE_OBSERVATIONS'] = \
#     os.path.join(home_dir, 'cascasde_test_hst_data')
# os.environ['CASCADE_INITIALIZATION_FILE_PATH'] = \
#     os.path.join(home_dir, 'cascade_test_ini_files')

# from cascade.initialize import cascade_default_data_path
# from cascade.initialize import cascade_default_initialization_path
# from cascade.exoplanet_tools import parse_database
# from cascade.exoplanet_tools import extract_exoplanet_data

try:
    import colorama
    colorama.init()
except ImportError:
    colorama = None

try:
    from termcolor import colored
except ImportError:
    colored = None

__all__ = ['built_local_hst_archive']


def log(string, color, font="slant", figlet=False):
    if colored:
        if not figlet:
            six.print_(colored(string, color))
        else:
            six.print_(colored(figlet_format(
                string, font=font), color))
    else:
        six.print_(string)


# ####################### START ###########################

@click.command()
@click.option('--init_path',
              '-ip',
              nargs=1,
              type=click.STRING,
              help='Path to the directory containing the configuration files.'
                   'If not specified, the value set by the environment '
                   'variable CASCADE_INITIALIZATION_FILE_PATH is used or if '
                   'neither is set it defaults to the CASCADe default value '
                   'of the CASCADe distribution. This path setting can be '
                   'relative to the absolute path given by the '
                   'CASCADE_INITIALIZATION_FILE_PATH environment variable.',
              )
@click.option('--data_path',
              '-dp',
              nargs=1,
              type=click.STRING,
              help='Path to the data (observations, calibration files, '
                   'exoplanet catalogs, archive databases) needed by CASCADe '
                   'to function. If not set, the value set by the environment '
                   'variable CASCADE_DATA_PATH is used or if neither is set '
                   'it defaults to the CASCADe default value of the CASCADe '
                   'distribution.',
              )
@click.option('--no_warnings',
              '-nw',
              is_flag=True,
              help='If set no warning messages are printed to stdev. '
                   'Default is True',
              )
def built_local_hst_archive(init_path, data_path, no_warnings):
    if no_warnings:
        os.environ["CASCADE_WARNINGS"] = 'off'
    else:
        os.environ["CASCADE_WARNINGS"] = 'on'
    if init_path is not None:
        if os.path.isabs(init_path):
            os.environ["CASCADE_INITIALIZATION_FILE_PATH"] = init_path
        else:
            try:
                base_ini_path = os.environ["CASCADE_INITIALIZATION_FILE_PATH"]
                os.environ["CASCADE_INITIALIZATION_FILE_PATH"] = \
                    os.path.join(base_ini_path, init_path)
            except KeyError:
                log("Relative init_path given without setting the "
                    "CASCADE_INITIALIZATION_FILE_PATH environment variable "
                    "Stooping script", "red")
                sys.exit()
    if data_path is not None:
        os.environ["CASCADE_DATA_PATH"] = data_path
#    import cascade
    from cascade.initialize import cascade_default_data_path
    from cascade.initialize import cascade_default_initialization_path
    from build_archive import return_exoplanet_catalogs
    from build_archive import read_config_file
    from build_archive import return_all_hst_planets
    from build_archive import return_hst_data_calalog_keys
    from build_archive import fill_system_parameters
    from build_archive import create_configuration
    from build_archive import print_parser_content
    from build_archive import long_substr
    from build_archive import return_header_info
    from build_archive import fill_config_parameters
    from build_archive import save_observations

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


if __name__ == '__main__':
    built_local_hst_archive()
    
