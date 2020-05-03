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
import pickle
import click
import sys
import six
from pyfiglet import figlet_format
from click import Option, UsageError
import time


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


class MutuallyExclusiveOption(Option):
    def __init__(self, *args, **kwargs):
        self.mutually_exclusive = set(kwargs.pop('mutually_exclusive', []))
        help = kwargs.get('help', '')
        if self.mutually_exclusive:
            ex_str = ', '.join(self.mutually_exclusive)
            kwargs['help'] = help + (
                ' NOTE: This argument is mutually exclusive with '
                ' arguments: [' + ex_str + '].'
            )
        super(MutuallyExclusiveOption, self).__init__(*args, **kwargs)

    def handle_parse_result(self, ctx, opts, args):
        if self.mutually_exclusive.intersection(opts) and self.name in opts:
            raise UsageError(
                "Illegal usage: `{}` is mutually exclusive with "
                "arguments `{}`.".format(
                    self.name,
                    ', '.join(self.mutually_exclusive)
                )
            )

        return super(MutuallyExclusiveOption, self).handle_parse_result(
            ctx,
            opts,
            args
        )


def check_path_option(new_path, environent_variable, message):
    if new_path is not None:
        if os.path.isabs(new_path):
            os.environ[environent_variable] = new_path
        else:
            try:
                base_new_path = os.environ[environent_variable]
                os.environ[environent_variable] = \
                    os.path.join(base_new_path, new_path)
            except KeyError:
                log(message, "red")
                sys.exit()


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
@click.option('--save_path',
              '-sp',
              nargs=1,
              type=click.STRING,
              help='Path to the directory where CASCADe saves results '
                   '(plots, extracted planetary spectra). If not specified, '
                   'the value set by the environment variable '
                   'CASCADE_SAVE_PATH is used, or if neither is set it '
                   'defaults to the CASCADe default value of the CASCADe '
                   'distribution.',
              )
@click.option('--no_warnings',
              '-nw',
              is_flag=True,
              default=False,
              help='If set no warning messages are printed to stdev. '
                   'Default is False',
              )
@click.option('--primary_exoplanet_catalog',
              '-pc',
              nargs=1,
              type=click.STRING,
              default='NASAEXOPLANETARCHIVE',
              help='The name of the primary exoplanet catalog used to create '
                   'the object.ini file'
              )
@click.option('--list_all_planets',
              '-lap', cls=MutuallyExclusiveOption,
              is_flag=True,
              default=False,
              help='Return all exoplanet names in Archive file.',
              mutually_exclusive=['list_catalog_id', 'visits',
                                  'all_visits_planet', 'download_all_data'],
              )
@click.option('--list_catalog_id',
              '-lci', cls=MutuallyExclusiveOption,
              nargs=1,
              type=click.STRING,
              help='List the observation calatog id of the planet which '
                   'name is given as an argument',
              mutually_exclusive=['list_all_planets', 'visits',
                                  'all_visits_planet', 'download_all_data']
              )
@click.option('--visits',
              '-v', cls=MutuallyExclusiveOption,
              multiple=True,
              help="List observation calatog id's for which the "
                   "initialization files will be created and the archive "
                   "data will be downloaded",
              mutually_exclusive=['list_all_planets', 'list_catalog_id',
                                  'all_visits_planet', 'download_all_data']
              )
@click.option('--all_visits_planet',
              '-avp', cls=MutuallyExclusiveOption,
              nargs=1,
              type=click.STRING,
              help="List observation calatog id's for a certain planet for "
                   "which the initialization files will be created and the "
                   "archive data will be downloaded",
              mutually_exclusive=['list_all_planets', 'list_catalog_id',
                                  'visits', 'download_all_data']
              )
@click.option('--download_all_data',
              '-dad', cls=MutuallyExclusiveOption,
              is_flag=True,
              default=False,
              help='Download all the archived observations to local disk.',
              mutually_exclusive=['list_all_planets', 'list_catalog_id',
                                  'visits', 'all_visits_planet'],
              )
@click.option('--create_ini_files_only',
              '-cio',
              is_flag=True,
              help='If set only the ini files are created. '
                   'Default is False'
              )
@click.option('--skip_existing_data',
              '-sed',
              is_flag=True,
              default=False,
              help='If set already downloaded data is skipped. '
                   'Default is False'
              )
@click.option('--skip_existing_initalization_files',
              '-sei',
              is_flag=True,
              default=False,
              help='If set already created initialization files are skipped. '
                   'Default is False'
              )
def built_local_hst_archive(init_path, data_path, save_path, no_warnings,
                            primary_exoplanet_catalog, list_all_planets,
                            list_catalog_id, visits, all_visits_planet,
                            download_all_data, create_ini_files_only,
                            skip_existing_data,
                            skip_existing_initalization_files):
    """
    Build local HST archive.

    This function creates a local copy of exoplanet observations in the
    MAST archive.
    """
    if no_warnings:
        os.environ["CASCADE_WARNINGS"] = 'off'
    else:
        os.environ["CASCADE_WARNINGS"] = 'on'

    if init_path is not None:
        check_path_option(init_path, "CASCADE_INITIALIZATION_FILE_PATH",
                          "Relative init_path given without setting the "
                          "CASCADE_INITIALIZATION_FILE_PATH environment "
                          "variable Stopping script")

    if data_path is not None:
        check_path_option(data_path, "CASCADE_DATA_PATH",
                          "Relative data_path given without setting the "
                          "CASCADE_DATA_PATH environment "
                          "variable Stopping script")

    if save_path is not None:
        check_path_option(save_path, "CASCADE_SAVE_PATH",
                          "Relative save_path given without setting the "
                          "CASCADE_SAVE_PATH environment "
                          "variable Stopping script")

    import cascade
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
    from build_archive import IniFileParser
    from build_archive import create_bash_script

    log('CASCADe', color="blue", figlet=True)
    log("version {}, Copyright (C) 2020 "
        "EXOPANETS_A H2020 Program".format(cascade.__version__), "blue")
    log("The initialization file directory is set to: "
        "{}".format(cascade_default_initialization_path), "green")
    log("The data directory is set to: "
        "{}".format(cascade_default_data_path), "green")

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
    if list_all_planets:
        log("The observations catalog contains data for the following "
            "planets: \n {}".format(return_all_hst_planets(hst_data_catalog)),
            'green')
        sys.exit()
    if list_catalog_id is not None:
        log("For planet {}, the following observations catalog keys have been "
            "found: {}".format(list_catalog_id,
                               return_hst_data_calalog_keys(list_catalog_id,
                                                            hst_data_catalog)),
            'green')
        sys.exit()

    if all_visits_planet is not None:
        visits = return_hst_data_calalog_keys(all_visits_planet,
                                              hst_data_catalog)
        log("For planet {}, the following observations catalog keys have been "
            "found: {}".format(all_visits_planet, visits),
            'green')

    if (visits is None) & (all_visits_planet is None) & \
            (not download_all_data):
        log("Warning, neither visits and all_visits_planet are defined "
            "either speciefy on of the other of use the --all_data option "
            "to download the entire archive.", "red")
        sys.exit()

    if download_all_data:
        if click.confirm("Selected --download_all_data option. Do you want "
                         "to continue downloading the entier archive?"):
            visits = hst_data_catalog.keys()
        else:
            log("Aborting downloading the entier archive", "red")
            sys.exit()

    if len(visits) == 0:
        log("Warning, no visits found, check search or path settings", "red")
        sys.exit()

    # explanet data catalogs
    catalogs_dictionary = return_exoplanet_catalogs()

    # ########### HERE LOOP STARTS #####################
    for visit in visits:
        PLANET_NAME = hst_data_catalog[visit]['planet']
        log("Target: {}, Visit: {}".format(PLANET_NAME, visit), 'green')

        # ############  some logic  ##############
        # ## Define the observation typpe and skips loop if problem ##
        if hst_data_catalog[visit]['observation'] == 'transit':
            OBS_TYPE = 'TRANSIT'
        elif hst_data_catalog[visit]['observation'] == 'eclipse':
            OBS_TYPE = 'ECLIPSE'
        else:
            OBS_TYPE = 'PROBLEM'
        if OBS_TYPE == 'PROBLEM':
            log('Skipping : {}'.format(hst_data_catalog[visit]['observation']),
                'red')
            continue
        # get the data files from hst databes
        data_files = hst_data_catalog[visit]['observations_id_ima'].split(',')
        cal_data_files = hst_data_catalog[visit]['calibrations_id'].split(',')
        data_file_id = [file.split('_')[0] for file in data_files]
        cal_data_file_id = [file.split('_')[0] for file in cal_data_files]
        # get or calculate all parameters needed for the observations and
        # instrument sections in the .ini file
        OBS_ID = long_substr(data_file_id)

        # ################ CREATE OBJECT.INI ##############################
        # read object.conf file
        object_config_dict = \
            read_config_file(OBJECT_CONFIGURATION_FILE, TEMPLATES_DIR)
        # fill with exoplanet catalog data
        object_config_dict = \
            fill_system_parameters(PLANET_NAME, catalogs_dictionary,
                                   object_config_dict,
                                   primary_exoplanet_catalog)
        # read catalog.conf file and fill parameters with values.
        # set which catalog to use
        catalog_name = primary_exoplanet_catalog
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
        # fill dictionary with parameters to be filled into templates
        create_timeseries_namespace_dict['cascade_save_path'] = \
            PLANET_NAME+'_'+OBS_ID
        create_timeseries_namespace_dict['observations_type'] = OBS_TYPE
        create_timeseries_namespace_dict['observations_target_name'] = \
            PLANET_NAME+'_'+OBS_ID
        create_timeseries_namespace_dict['observations_id'] = OBS_ID
        create_timeseries_namespace_dict.update(
            **return_header_info(data_files[0], cal_data_files[0])
            )
        # create list with the configuration files to use
        configuration_file_list = [CASCADE_CONFIGURATION_FILE,
                                   PROCESSING_CONFIGURATION_FILE,
                                   INSTRUMENT_CONFIGURATION_FILE,
                                   OBSERVATIONS_CONFIGURATION_FILE]
        # create initialization file parser object
        IFP = \
            IniFileParser(configuration_file_list,
                          EXTRACT_TIMESERIES_CONFIGURATION_TEMPLATE,
                          create_timeseries_namespace_dict,
                          TEMPLATES_DIR)
        extract_timeseries_parser = IFP.return_parser()
        print_parser_content(extract_timeseries_parser)

        # ############### CREATE CALIBRATE_PLANET_SPECTRUM.INI ################
        # update parameters for timeseries of 1D spectra
        cal_planet_spec_namespace_dict = \
            create_timeseries_namespace_dict.copy()
        cal_planet_spec_namespace_dict['processing_nextraction'] = 1
        cal_planet_spec_namespace_dict['observations_data_product'] = 'COE'
        cal_planet_spec_namespace_dict['observations_has_background'] = False
        cal_planet_spec_namespace_dict['observations_mode'] = 'STARING'
        cal_planet_spec_namespace_dict['observations_data'] = 'SPECTRUM'

        # update which configureation files to use
        configuration_file_list = [CASCADE_CONFIGURATION_FILE,
                                   PROCESSING_CONFIGURATION_FILE,
                                   CPM_CONFIGURATION_FILE,
                                   MODEL_CONFIGURATION_FILE,
                                   INSTRUMENT_CONFIGURATION_FILE,
                                   OBSERVATIONS_CONFIGURATION_FILE]
        # create initialization file parser object
        IFP = \
            IniFileParser(configuration_file_list,
                          CALIBRATE_PLANET_SPECTRUM_CONFIGURATION_TEMPLATE,
                          cal_planet_spec_namespace_dict,
                          TEMPLATES_DIR)
        calibrate_planet_spectrum_parser = IFP.return_parser()
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
        if not (skip_existing_initalization_files & os.path.isfile(
            os.path.join(initialization_save_path,
                         configuration_object_filename))):
            with open(os.path.join(initialization_save_path,
                                   configuration_object_filename),
                      'w') as configfile:
                object_parser.write(configfile)

        configuration_extract_timeseries_filename = \
            configuration_base_filename+"_extract_timeseries.ini"

        if not (skip_existing_initalization_files & os.path.isfile(
            os.path.join(initialization_save_path,
                         configuration_extract_timeseries_filename))):
            with open(os.path.join(
                    initialization_save_path,
                    configuration_extract_timeseries_filename),
                    'w') as configfile:
                extract_timeseries_parser.write(configfile)

        configuration_calibrate_planet_spectrum_filename = \
            configuration_base_filename+"_calibrate_planet_spectrum.ini"
        if not (skip_existing_initalization_files & os.path.isfile(
            os.path.join(initialization_save_path,
                         configuration_calibrate_planet_spectrum_filename))):
            with open(os.path.join(
                    initialization_save_path,
                    configuration_calibrate_planet_spectrum_filename),
                    'w') as configfile:
                calibrate_planet_spectrum_parser.write(configfile)

        # # ################# GET ARCHIVE DATA ######################
        if not create_ini_files_only:
            save_observations(data_files, cal_data_files,
                              extract_timeseries_parser,
                              skip_existing=skip_existing_data)
        else:
            log("create_ini_files_only flag is True, not downloading data",
                "red")

        # ###################### crete bash script #######################
        create_bash_script(visit, extract_timeseries_parser)


if __name__ == '__main__':
    start_time = time.time()
    built_local_hst_archive()
    elapsed_time = time.time() - start_time
    log('elapsed time: {}'.format(elapsed_time), "green")
